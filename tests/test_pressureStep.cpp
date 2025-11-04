#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath> // For std::sin, std::cos, std::pow

// --- Include headers ---
#include "simulation/pressureStep.hpp" // Object under test
#include "simulation/SimulationContext.hpp"
#include "core/Fields.hpp"
#include "numerics/derivatives.hpp" // For comparison


/**
 * @brief Test Fixture for unit and correctness tests for PressureStep.
 * Uses a 10x10x10 grid with dx=1.0 and non-zero initial conditions.
 */
class PressureStepTest : public ::testing::Test {
protected:
    // Test objects
     const double dt = 0.1;
    std::shared_ptr<Grid> grid;
    SimulationContext context_;
    std::unique_ptr<PressureStep> pressureStep;


    // --- CONSTRUCTOR ---
    PressureStepTest()
        : // 1. Create the grid (10x10x10, dx=1.0)
          grid(std::make_shared<Grid>(10, 10, 10, 1.0, 1.0, 1.0)),
          // 2. Initialize the context_
          context_{ grid, {dt, 1.0} } // dt=0.1, totalTime=1.0
    {
        // Constructor body is empty. SetUp() runs next.
    }

    // --- SETUP ---
    // Runs before each test in this fixture
    void SetUp() override {
        
        std::vector<Field::Scalar> zeros(grid->size(), 0.0);
        std::vector<Field::Scalar> u_x(grid->size()), u_y(grid->size()), u_z(grid->size());
        
        // --- Set Initial Conditions ---
        
        // 1. Pressure: p(x,y,z) = 0.0
        //    PressureStep will add 'pcr' to this.
        std::vector<Field::Scalar> p_field(grid->size(), 0.0);

        // 2. Velocity: u(x,y,z) = [sin(x), cos(y), sin(z)]
        //    Whose divergence is: cos(x) - sin(y) + cos(z)
        for (size_t k = 0; k < grid->Nz; ++k) {
            for (size_t j = 0; j < grid->Ny; ++j) {
                for (size_t i = 0; i < grid->Nx; ++i) {
                    size_t index = i + grid->Nx * (j + grid->Ny * k);
                    double x = static_cast<double>(i) * grid->dx;
                    double y = static_cast<double>(j) * grid->dy;
                    double z = static_cast<double>(k) * grid->dz;

                    u_x[index] = std::sin(x);
                    u_y[index] = std::cos(y);
                    u_z[index] = std::sin(z);
                }
            }
        }
        
        // Setup fields in state
        context_.state.p.setup(grid, p_field);
        context_.state.u.setup(grid, u_x, u_y, u_z); // Current 'u'
        
        // Remove unneeded fields (eta, zeta, uOld, etc.)

        // 4. Create the object under test
        pressureStep = std::make_unique<PressureStep>(context_);
        
        // NOTE: Your initializeWorkspaceFields only inits psi and divU.
        // If phi and pcr are also members, you might need to init them here.
        // pressureStep->phi.setup(grid, zeros);
        // pressureStep->pcr.setup(grid, zeros);
    }

    // --- HELPER METHODS ---
protected:
    double getDivU(size_t i, size_t j, size_t k) {
        // Assuming 'divU' is a member (public or friend-accessible)
        return pressureStep->divU(i, j, k);
    }
    
    double getPressure(size_t i, size_t j, size_t k) {
        return context_.state.p(i, j, k);
    }
};


// --- The Tests ---

TEST_F(PressureStepTest, ConstructorAndInit) {
    ASSERT_NE(pressureStep, nullptr);
}

TEST_F(PressureStepTest, RunDoesNotCrash) {
    // This test verifies that the run() function completes without errors
    ASSERT_NO_THROW(pressureStep->run());
}

TEST_F(PressureStepTest, Run_CalculatesDivU_Correctness) {
    // This test verifies that the divergence calculated inside run()
    // is correct.
    
    // --- Execution ---
    ASSERT_NO_THROW(pressureStep->run());

    // --- Verify Output (at internal node 5,5,5) ---
    const size_t i = 5, j = 5, k = 5;
    const double x = i*grid->dx, y = j*grid->dy, z = k*grid->dz;
    const double dx = grid->dx, dy = grid->dy, dz = grid->dz;

    // u.x = sin(x) -> d/dx = (sin(x+dx) - sin(x-dx)) / (2*dx)
    double du_dx = (std::sin(x+dx) - std::sin(x-dx)) / (2.0*dx);
    // u.y = cos(y) -> d/dy = (cos(y+dy) - cos(y-dy)) / (2*dy)
    double du_dy = (std::cos(y+dy) - std::cos(y-dy)) / (2.0*dy);
    // u.z = sin(z) -> d/dz = (sin(z+dz) - sin(z-dz)) / (2*dz)
    double du_dz = (std::sin(z+dz) - std::sin(z-dz)) / (2.0*dz);
    
    // div(u) = du_dx + du_dy + du_dz
    double expected_divU = du_dx + du_dy + du_dz;
    
    double divU_computed = getDivU(i, j, k);
    EXPECT_NEAR(divU_computed, expected_divU, 1e-9);
}

TEST_F(PressureStepTest, Run_WithZeroDivergence_PressureIsUnchanged) {
    // This test sets u = [0,0,0], which has zero divergence.
    // The resulting pcr should be 0, and the final pressure
    // should be equal to the initial pressure.

    // --- Setup (Override) ---
    std::vector<Field::Scalar> zeros(grid->size(), 0.0);
    std::vector<Field::Scalar> ones(grid->size(), 1.0);
    
    // Set P_initial = 1.0
    context_.state.p.setup(grid, ones);
    // Set U = 0
    context_.state.u.setup(grid, zeros, zeros, zeros);

    // Recreate the test object with the new context_
    pressureStep = std::make_unique<PressureStep>(context_);
    // ...potential init of phi/pcr...
    // pressureStep->phi.setup(grid, zeros);
    // pressureStep->pcr.setup(grid, zeros);


    // --- Execution ---
    pressureStep->run();

    // --- Verification ---
    // p_final = p_initial + pcr
    // If div(u) = 0, then pcr should be 0.
    // p_final should be 1.0
    
    double p_final = getPressure(5, 5, 5); // Check one point
    EXPECT_NEAR(p_final, 1.0, 1e-9);
    
    // Check that divU is actually 0
    double divU_computed = getDivU(5, 5, 5);
    EXPECT_NEAR(divU_computed, 0.0, 1e-9);
}


/**
 * @brief Test Fixture for testing solver robustness with polynomial fields.
 */
class PressureStepRobustnessTest : public ::testing::Test {
protected:
    const size_t N = 10;
    std::shared_ptr<Grid> grid;
    SimulationContext context_;
    std::unique_ptr<PressureStep> pressureStep;

    PressureStepRobustnessTest()
        : // 1. Grid 10x10x10, dx=dy=dz=0.1
          grid(std::make_shared<Grid>(N, N, N, 0.1, 0.1, 0.1)),
          // 2. Context: dt=0.01
          context_{ grid, {0.01, 1.0} }
    {}

   void SetUp() override {
        
        const size_t totalSize = N * N * N;
        std::vector<Field::Scalar> p_field(totalSize);
        std::vector<Field::Scalar> u_x(totalSize), u_y(totalSize), u_z(totalSize);

        // Create non-zero polynomial fields
        for (size_t k = 0; k < N; ++k) {
            for (size_t j = 0; j < N; ++j) {
                for (size_t i = 0; i < N; ++i) {
                    size_t index = i + N * (j + N * k);
                    double x = static_cast<double>(i) * grid->dx;
                    double y = static_cast<double>(j) * grid->dy;
                    double z = static_cast<double>(k) * grid->dz;

                    // p = x*y + z
                    p_field[index] = x*y + z;
                    
                    // u = [x^2, y^2, z^2] -> div(u) = 2x + 2y + 2z
                    u_x[index] = x*x;
                    u_y[index] = y*y;
                    u_z[index] = z*z;
                }
            }
        }

        // --- Set up the context_ ---
        context_.state.p.setup(grid, p_field);
        context_.state.u.setup(grid, u_x, u_y, u_z); // Initial u
        
        pressureStep = std::make_unique<PressureStep>(context_);
    }

    // Helper to check for NaN or Inf in a scalar field
    void checkFieldFinite(const Field& field, const std::string& fieldName) {
        // Assuming Field::Scalar has an operator()
        
        for (size_t k = 0; k < grid->Nz; ++k) {
            for (size_t j = 0; j < grid->Ny; ++j) {
                for (size_t i = 0; i < grid->Nx; ++i) {
                     ASSERT_TRUE(std::isfinite(field(i, j, k)))
                        << "Found NaN or Inf in " << fieldName 
                        << " at index (" << i << "," << j << "," << k << ")";
                }
            }
        }
    }
};

/**
 * @brief Robustness test for the full PressureStep::run()
 *
 * This test uses complex, non-zero polynomial initial conditions
 * and checks that the solver produces valid, finite numbers
 * (i.e., no NaN/Inf) in the final pressure field.
 */
TEST_F(PressureStepRobustnessTest, Run_WithPolynomialFields_ProducesFiniteOutput)
{
    // --- Execution ---
    ASSERT_NO_THROW(pressureStep->run());

    // --- Verification ---
    // We only check that the final pressure is stable
    checkFieldFinite(context_.state.p, "context_.state.p");
}