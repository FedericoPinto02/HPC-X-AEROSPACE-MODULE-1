#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>
#include <functional>

#include "simulation/viscousStep.hpp"
#include "simulation/SimulationContext.hpp"
#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"

class ViscousStepTest : public ::testing::Test {
protected:
    const double dt = 0.1;
    const double nu = 1.0;

    // Grid managed via shared_ptr to match Field's internal storage requirement
    std::shared_ptr<Grid> grid;

    // Context
    SimulationData data_;
    ParallelizationSettings parallel_;

    std::unique_ptr<ViscousStep> viscousStep;

    ViscousStepTest() {
        // Initialize 10x10x10 grid with spacing 1.0
        grid = std::make_shared<Grid>(10, 10, 10, 1.0, 1.0, 1.0);

        // Setup SimulationData context
        data_.grid = grid;
        data_.dt = dt;
        data_.nu = nu;
        data_.currTime = 0.0;

        parallel_.schurDomains = 1;
    }

    void SetUp() override {
        // Define constant Permeability k = 1.0
        data_.inv_k.setup(grid, [](double, double, double, double) { return 1.0; });
        data_.inv_k.populate(0.0);

        // Define External Force f = [0,0,0]
        data_.f.setup(grid, ZERO_FUNC, ZERO_FUNC, ZERO_FUNC);
        data_.f.populate(0.0);

        // Define Pressure: p(x,y,z) = x + 2y + 3z
        auto funcP = [](double x, double y, double z, double) {
            return x + 2.0 * y + 3.0 * z;
        };
        data_.p.setup(grid, funcP);
        data_.p.populate(0.0);

        // Define Pressure predictor: same as p
        data_.predictor.setup(grid, funcP);
        data_.predictor.populate(0.0);

        // Define Velocity: u = [sin(x), sin(y), sin(z)]
        auto sinX = [](double x, double, double, double) { return std::sin(x); };
        auto sinY = [](double, double y, double, double) { return std::sin(y); };
        auto sinZ = [](double, double, double z, double) { return std::sin(z); };
        data_.u.setup(grid, sinX, sinY, sinZ);
        data_.u.populate(0.0);

        // Define Eta (intermediate velocity): eta = [cos(x), cos(y), cos(z)]
        auto cosX = [](double x, double, double, double) { return std::cos(x); };
        auto cosY = [](double, double y, double, double) { return std::cos(y); };
        auto cosZ = [](double, double, double z, double) { return std::cos(z); };
        data_.eta.setup(grid, cosX, cosY, cosZ);
        data_.eta.populate(0.0);

        // Define Zeta (intermediate velocity): zeta = [x^2, y^2, z^2]
        auto sqX = [](double x, double, double, double) { return x * x; };
        auto sqY = [](double, double y, double, double) { return y * y; };
        auto sqZ = [](double, double, double z, double) { return z * z; };
        data_.zeta.setup(grid, sqX, sqY, sqZ);
        data_.zeta.populate(0.0);

        // Define Analytical Boundary Functions
        // Required for LinearSys solver to evaluate BCs at specific coordinates
        data_.bcu = [](double x, double, double, double) { return std::sin(x); };
        data_.bcv = [](double, double y, double, double) { return std::sin(y); };
        data_.bcw = [](double, double, double z, double) { return std::sin(z); };

        viscousStep = std::make_unique<ViscousStep>(data_, parallel_);
    }

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

    void setGConstant(double val) {
        auto funcVal = [val](double, double, double, double) { return val; };
        viscousStep->g.setup(grid, funcVal, funcVal, funcVal);
        viscousStep->g.populate(0.0);
    }
};

TEST_F(ViscousStepTest, ConstructorAndInit) {
    ASSERT_NE(viscousStep, nullptr);
}

TEST_F(ViscousStepTest, RunDoesNotCrash) {
// Triggers LinearSys::fillSystemVelocity which calls analytical BCs
    ASSERT_NO_THROW(viscousStep->run());
}

TEST_F(ViscousStepTest, ComputeG_Correctness) {
    const size_t i = 5, j = 5, k_ = 5;
    const double dx = grid->dx;
    const double k_permeability = 1.0;

    ASSERT_NO_THROW(callComputeG());

    // The G term for the X-component lives at (i+0.5)
    double x = grid->to_x(i, GridStaggering::FACE_CENTERED, Axis::X);
    double x_prev = grid->to_x(i - 1, GridStaggering::FACE_CENTERED, Axis::X);
    double x_next = grid->to_x(i + 1, GridStaggering::FACE_CENTERED, Axis::X);

    double f_term = 0.0;
    double gradP_term = 1.0; // d/dx(x + 2y + 3z) = 1
    double u_term = (nu / k_permeability) * std::sin(x);

// Laplacian approximation
    double dxx_term = (std::cos(x_next) - 2.0 * std::cos(x) + std::cos(x_prev)) / (dx * dx);
    double dyy_term = 0.0; // Assume constant in Y for this test
    double dzz_term = 0.0; // Assume constant in Z for this test
    double deriv_term = nu * (dxx_term + dyy_term + dzz_term);

// Expected G formula from implementation
    double expected_g_x = f_term - gradP_term - u_term + deriv_term;

    double g_x_computed = getG(Axis::X, i, j, k_);
    EXPECT_NEAR(g_x_computed, expected_g_x, 1e-9);
}

TEST_F(ViscousStepTest, ComputeXi_Correctness) {
    const size_t i = 5, j = 5, k_ = 5;
    const double x = grid->to_x(i, GridStaggering::FACE_CENTERED, Axis::X);
    double u_val = std::sin(x);

    setGConstant(10.0);

    callComputeXi();

    double k_val = 1.0;
    double beta = 1.0 + (dt * nu * 0.5 / k_val);

    double expected_xi_x = u_val + (dt / beta) * 10.0;

    double xi_x_computed = getXi(Axis::X, i, j, k_);
    EXPECT_NEAR(xi_x_computed, expected_xi_x, 1e-9);
}

class ViscousStepRobustnessTest : public ::testing::Test {
protected:
    const size_t N = 10;
    std::shared_ptr<Grid> grid;
    SimulationData data_;
    ParallelizationSettings parallel_;
    std::unique_ptr<ViscousStep> viscousStep;

    ViscousStepRobustnessTest() {
        grid = std::make_shared<Grid>(N, N, N, 0.1, 0.1, 0.1);
        data_.grid = grid;
        data_.dt = 0.1;
        data_.nu = 1.0e-3;
        data_.currTime = 0.0;
        parallel_.schurDomains = 1;
    }

    void SetUp() override {
        // Define Permeability k = 1.0
        data_.inv_k.setup(grid, [](double, double, double, double) { return 1.0; });
        data_.inv_k.populate(0.0);

        // Define Force f = 10.0
        auto funcF = [](double, double, double, double) { return 10.0; };
        data_.f.setup(grid, funcF, funcF, funcF);
        data_.f.populate(0.0);

        // Define Pressure p = x*y + z
        data_.p.setup(grid, [](double x, double y, double z, double) { return x * y + z; });
        data_.p.populate(0.0);

        // Define Pressure predictor: same as p
        data_.predictor.setup(grid, [](double x, double y, double z, double) { return x * y + z; });
        data_.predictor.populate(0.0);


        // Define Velocity u = [x^2, y^2, z^2]
        auto sqX = [](double x, double, double, double) { return x * x; };
        auto sqY = [](double, double y, double, double) { return y * y; };
        auto sqZ = [](double, double, double z, double) { return z * z; };
        data_.u.setup(grid, sqX, sqY, sqZ);
        data_.u.populate(0.0);

        // Define Eta = [y^2, z^2, x^2]
        data_.eta.setup(grid, sqY, sqZ, sqX);
        data_.eta.populate(0.0);

        // Define Zeta = [z^2, x^2, y^2]
        data_.zeta.setup(grid, sqZ, sqX, sqY);
        data_.zeta.populate(0.0);

        // Define Analytical Boundary Functions
        data_.bcu = sqX;
        data_.bcv = sqY;
        data_.bcw = sqZ;

        viscousStep = std::make_unique<ViscousStep>(data_, parallel_);
    }

    void checkFieldFinite(const Field &field, const std::string &fieldName) {
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

TEST_F(ViscousStepRobustnessTest, Run_WithPolynomialFields_ProducesFiniteOutput) {
    ASSERT_NO_THROW(viscousStep->run());

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