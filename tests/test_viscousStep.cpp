#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath> // For std::sin, std::cos, std::pow

// --- Include headers ---
#include "simulation/viscousStep.hpp"
#include "simulation/SimulationContext.hpp"
#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"

// Helper function to fill a vector with a sin(a*coord)
void fillSines(std::vector<Field::Scalar>& x, 
               std::vector<Field::Scalar>& y, 
               std::vector<Field::Scalar>& z, 
               std::shared_ptr<const Grid> grid)
{
    for (size_t k = 0; k < grid->Nz; ++k) {
        for (size_t j = 0; j < grid->Ny; ++j) {
            for (size_t i = 0; i < grid->Nx; ++i) {
                size_t index = i + grid->Nx * (j + grid->Ny * k);
                double x_coord = static_cast<double>(i) * grid->dx;
                double y_coord = static_cast<double>(j) * grid->dy;
                double z_coord = static_cast<double>(k) * grid->dz;
                x[index] = std::sin(x_coord);
                y[index] = std::sin(y_coord);
                z[index] = std::sin(z_coord);
            }
        }
    }
}


/**
 * @brief Test Fixture for unit and correctness tests.
 * Uses a 10x10x10 grid with dx=1.0 and non-zero initial conditions.
 */
class ViscousStepTest : public ::testing::Test {
protected:
    // Test objects

    const double dt = 0.1;
    std::shared_ptr<Grid> grid;
    SimulationData data_;
    std::unique_ptr<ViscousStep> viscousStep;

    const double nu = 1.0;

    // --- CONSTRUCTOR ---
    ViscousStepTest()
        : // 1. Create the grid (10x10x10, dx=1.0)
          grid(std::make_shared<Grid>(10, 10, 10, 1.0, 1.0, 1.0)),
          // 2. Initialize the data_
          data_{ grid } // dt=0.1, totalTime=1.0
    {
        // Constructor body is empty. SetUp() runs next.
    }

    // --- SETUP ---
    // Runs before each test in this fixture
    void SetUp() override {
        
        std::vector<Field::Scalar> zeros(grid->size(), 0.0);
        std::vector<Field::Scalar> ones(grid->size(), 1.0); 

        // Set constants
        data_.nu = nu;
        data_.dt = dt;
        data_.k.setup(grid, ones); // k = 1.0 everywhere
        data_.f.setup(grid, zeros, zeros, zeros); // f = [0,0,0]

        // --- Set "Reasonable" Non-Zero Initial Conditions ---
        
        // 1. Pressure: p(x,y,z) = x + 2y + 3z
        //    (Gradient should be [1, 2, 3])
        std::vector<Field::Scalar> p_field(grid->size());
        // 2. Velocity: uOld(x,y,z) = [sin(x), sin(y), sin(z)]
        std::vector<Field::Scalar> u_x(grid->size()), u_y(grid->size()), u_z(grid->size());
        // 3. Eta: etaOld(x,y,z) = [cos(x), cos(y), cos(z)]
        std::vector<Field::Scalar> eta_x(grid->size()), eta_y(grid->size()), eta_z(grid->size());
        // 4. Zeta: zetaOld(x,y,z) = [x^2, y^2, z^2]
        std::vector<Field::Scalar> zeta_x(grid->size()), zeta_y(grid->size()), zeta_z(grid->size());

        for (size_t k = 0; k < grid->Nz; ++k) {
            for (size_t j = 0; j < grid->Ny; ++j) {
                for (size_t i = 0; i < grid->Nx; ++i) {
                    size_t index = i + grid->Nx * (j + grid->Ny * k);
                    double x = static_cast<double>(i) * grid->dx;
                    double y = static_cast<double>(j) * grid->dy;
                    double z = static_cast<double>(k) * grid->dz;

                    p_field[index] = x + 2.0*y + 3.0*z;
                    
                    u_x[index] = std::sin(x);
                    u_y[index] = std::sin(y);
                    u_z[index] = std::sin(z);

                    eta_x[index] = std::cos(x);
                    eta_y[index] = std::cos(y);
                    eta_z[index] = std::cos(z);

                    zeta_x[index] = x*x;
                    zeta_y[index] = y*y;
                    zeta_z[index] = z*z;
                }
            }
        }
        
        // Setup fields
        data_.p.setup(grid, p_field);
        data_.u.setup(grid, u_x, u_y, u_z);
        data_.eta.setup(grid, eta_x, eta_y, eta_z);
        data_.zeta.setup(grid, zeta_x, zeta_y, zeta_z);

        // Initialize boundary condition fields
        data_.uBoundNew.setup(grid, u_x, u_y, u_z);
        data_.uBoundOld.setup(grid, u_x, u_y, u_z);
        // --- End Non-Zero ICs ---

        // 4. Create the object under test
        viscousStep = std::make_unique<ViscousStep>(data_);
    }

    // --- HELPER METHODS ---
protected:
    void callComputeG() {
        viscousStep->computeG();
    }
    
    void callComputeXi() {
        viscousStep->computeXi();
    }

    double getG(Axis axis, size_t i, size_t j, size_t k) {
        return viscousStep->g(axis, i, j, k);
    }

    double getXi(Axis axis, size_t i, size_t j, size_t k) {
        return viscousStep->xi(axis, i, j, k);
    }
    
    void setG(std::shared_ptr<const Grid> gridPtr, 
              std::vector<Field::Scalar>& x, 
              std::vector<Field::Scalar>& y, 
              std::vector<Field::Scalar>& z) {
        viscousStep->g.setup(gridPtr, x, y, z);
    }
};


// --- The Tests ---

TEST_F(ViscousStepTest, ConstructorAndInit) {
    ASSERT_NE(viscousStep, nullptr);
}

TEST_F(ViscousStepTest, RunDoesNotCrash) {
    // This test will only pass if:
    // 1. The Field::setup() bug (SIGSEGV) is fixed.
    // 2. The viscousStep.cpp NaN bug (local 'xi') is fixed.
    ASSERT_NO_THROW(viscousStep->run());
}

TEST_F(ViscousStepTest, ComputeG_Correctness) {
    // Formula: g = f - grad(p) - nu/k * u + nu * (dxx eta + dyy zeta + dzz u)
    // We check the 'X' component at an internal node (e.g., 5, 5, 5).
    const size_t i = 5, j = 5, k_ = 5;
    const double x = i*grid->dx, y = j*grid->dy, z = k_*grid->dz;

    // --- Execution ---
    // Must fix Field::setup bug before this will pass
    ASSERT_NO_THROW(callComputeG());

    // --- Verify Output (at internal node 5,5,5) ---
    
    // 1. f.x = 0 (from SetUp)
    double f_term = 0.0;
    
    // 2. grad(p).x = 1.0
    // (p[i+1] - p[i-1]) / (2*dx) = ((x+1)+2y+3z - ((x-1)+2y+3z)) / 2 = 2/2 = 1.0
    double gradP_term = 1.0;
    
    // 3. (nu/k)*u.x = (1.0/1.0) * sin(x)
    double u_term = (nu / 1.0) * std::sin(x);
    
    // 4. nu*(dxx(eta.x) + dyy(zeta.x) + dzz(u.x))
    // dxx(eta.x) = dxx(cos(x)) = (cos(x+1) - 2cos(x) + cos(x-1)) / 1^2
    double dxx_term = std::cos(x+1.0) - 2.0*std::cos(x) + std::cos(x-1.0);
    // dyy(zeta.x) = dyy(x^2) = 0
    double dyy_term = 0.0;
    // dzz(u.x) = dzz(sin(x)) = 0
    double dzz_term = 0.0;
    
    double deriv_term = nu * (dxx_term + dyy_term + dzz_term);
    
    // g.x = 0 - 1.0 - sin(x) + 1.0 * ( (cos(x+1)-2cos(x)+cos(x-1)) + 0 + 0 )
    double expected_g_x = f_term - gradP_term - u_term * 0.5 + deriv_term * 0.5;
    
    double g_x_computed = getG(Axis::X, i, j, k_);
    EXPECT_NEAR(g_x_computed, expected_g_x, 1e-9);
}

TEST_F(ViscousStepTest, ComputeXi_Correctness) {
    // Formula: xi = u + dt/beta * g
    //     beta = 1 + dt*nu / (2*k)
    
    // --- Setup Input ---
    // dt = 0.1, nu = 1.0, k = 1.0 (from SetUp)
    // uOld.x(5,5,5) = sin(5)
    const size_t i=5, j=5, k_=5;
    const double x = i*grid->dx;
    double u_val = std::sin(x);
    
    // 1. Manually set g.x = 10.0 everywhere
    std::vector<Field::Scalar> tens(grid->size(), 10.0);
    std::vector<Field::Scalar> zeros(grid->size(), 0.0);
    setG(grid, tens, zeros, zeros);

    // --- Execution ---
    callComputeXi(); // Call helper method

    // --- Verify Output ---
    // beta = 1 + (0.1 * 1.0) / (2 * 1.0) = 1.05
    double k_val = 1.0;
    double beta = 1.0 + (dt * nu * 0.5 / k_val);
    
    // xi.x = u.x + (dt / beta) * g.x
    // xi.x = sin(5.0) + (0.1 / 1.05) * 10.0
    double expected_xi_x = u_val + (dt / beta) * 10.0;
    
    double xi_x_computed = getXi(Axis::X, i, j, k_);
    EXPECT_NEAR(xi_x_computed, expected_xi_x, 1e-9);
}




/**
 * @brief Test Fixture for testing the closeViscousStep solver.
 * This inherits from ViscousStepTest to get the SetUp logic.
 */
class ViscousStepSolverTest : public ViscousStepTest {
protected:
    // Helper to call the private closeViscousStep method
    void callCloseViscousStep() {
        viscousStep->closeViscousStep();
    }
};


/**
 * @brief Test Fixture for testing solver robustness with polynomial fields.
 * This test uses complex, non-zero fields to ensure the solver
 * does not produce NaN or Inf values.
 */
class ViscousStepRobustnessTest : public ::testing::Test {
protected:
    const size_t N = 10;
    std::shared_ptr<Grid> grid;
    SimulationData data_;
    std::unique_ptr<ViscousStep> viscousStep;

   

    ViscousStepRobustnessTest()
        : // 1. Grid 5x5x5, dx=dy=dz=0.1
          grid(std::make_shared<Grid>(N, N, N, 0.1, 0.1, 0.1)),
          // 2. data_: dt=0.01
          data_{ grid }
    {}

   void SetUp() override {
        
        // Calcola la dimensione corretta manualmente
        const size_t totalSize = N * N * N;

        std::vector<Field::Scalar> ones(totalSize, 1.0);
        std::vector<Field::Scalar> f_field(totalSize, 10.0); // Forcing
        std::vector<Field::Scalar> p_field(totalSize);
        
        std::vector<Field::Scalar> u_x(totalSize), u_y(totalSize), u_z(totalSize);
        std::vector<Field::Scalar> eta_x(totalSize), eta_y(totalSize), eta_z(totalSize);
        std::vector<Field::Scalar> zeta_x(totalSize), zeta_y(totalSize), zeta_z(totalSize);


        // Create non-zero polynomial fields
        for (size_t k = 0; k < N; ++k) {
            for (size_t j = 0; j < N; ++j) {
                for (size_t i = 0; i < N; ++i) {
                    size_t index = i + N * (j + N * k);
                    // Use physical coordinates
                    double x = static_cast<double>(i) * grid->dx;
                    double y = static_cast<double>(j) * grid->dy;
                    double z = static_cast<double>(k) * grid->dz;

                    // p = x*y + z
                    p_field[index] = x*y + z;
                    
                    // uOld = [x^2, y^2, z^2]
                    u_x[index] = x*x;
                    u_y[index] = y*y;
                    u_z[index] = z*z;

                    // etaOld = [y^2, z^2, x^2]
                    eta_x[index] = y*y;
                    eta_y[index] = z*z;
                    eta_z[index] = x*x;
                    
                    // zetaOld = [z^2, x^2, y^2]
                    zeta_x[index] = z*z;
                    zeta_y[index] = x*x;
                    zeta_z[index] = y*y;
                }
            }
        }

        // --- Set up the data_ ---
        data_.dt = 0.1;
        data_.nu = 1.0e-3;
        data_.k.setup(grid, ones); // k = 1.0 (non-zero)
        data_.f.setup(grid, f_field, f_field, f_field); // f = [10, 10, 10]

        // Set state fields
        data_.p.setup(grid, p_field);
        data_.u.setup(grid, u_x, u_y, u_z); // Initial u
        data_.eta.setup(grid, eta_x, eta_y, eta_z); // Initial eta
        data_.zeta.setup(grid, zeta_x, zeta_y, zeta_z); // Initial zeta

        // Set non-zero boundary conditions
        data_.uBoundNew.setup(grid, u_x, u_y, u_z);
        data_.uBoundOld.setup(grid, u_x, u_y, u_z);
        
        viscousStep = std::make_unique<ViscousStep>(data_);
    }

    // Helper function to check a field for NaN or Inf
    void checkFieldFinite(const Field& field, const std::string& fieldName) {
        const auto& data = field.getData();
        for (size_t i = 0; i < data.size(); ++i) {
            ASSERT_TRUE(std::isfinite(data[i]))
                << "Found NaN or Inf in " << fieldName 
                << " at index " << i;
        }
    }
};

/**
 * @brief Robustness test for the full ViscousStep::run()
 *
 * This test uses complex, non-zero polynomial initial conditions
 * and checks that the solver produces valid, finite numbers
 * (i.e., no NaN or Inf) in the final output fields.
 *
 * This test will FAIL if:
 * 1. Your Field::setup() method does not set p_grid (causes SIGSEGV).
 * 2. Your viscousStep.cpp passes data_.xi to the solver (causes NaN).
 * 3. Your solver logic has a divide-by-zero or other instability.
 */
TEST_F(ViscousStepRobustnessTest, Run_WithPolynomialFields_ProducesFiniteOutput)
{
    // --- Execution ---
    // This will crash if the SIGSEGV bug in Field::setup is not fixed.
    // This will produce NaN if the 'xi' bug in viscousStep.cpp is not fixed.
    ASSERT_NO_THROW(viscousStep->run());

    // --- Verification ---
    // We don't check for specific values. We only check that the
    // simulation is stable and did not produce NaN or Inf.
    
    checkFieldFinite(data_.eta(Axis::X), "eta(X)");
    checkFieldFinite(data_.eta(Axis::Y), "eta(Y)");
    checkFieldFinite(data_.eta(Axis::Z), "eta(Z)");
    
    checkFieldFinite(data_.zeta(Axis::X), "zeta(X)");
    checkFieldFinite(data_.zeta(Axis::Y), "zeta(Y)");
    checkFieldFinite(data_.zeta(Axis::Z), "zeta(Z)");

    checkFieldFinite(data_.u(Axis::X), "u(X)");
    checkFieldFinite(data_.u(Axis::Y), "u(Y)");
    checkFieldFinite(data_.u(Axis::Z), "u(Z)");
}