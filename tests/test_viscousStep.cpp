#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>

#include "MpiEnvFixture.hpp"

#include "core/Fields.hpp"
#include "simulation/viscousStep.hpp"
#include "simulation/SimulationContext.hpp"

class ViscousStepTest : public ::testing::Test {
protected:
    const double dt = 0.1;
    const double nu = 1.0;

    // Grid managed via shared_ptr to match Field's internal storage requirement
    std::shared_ptr<const Grid> grid;

    // Context
    SimulationData data_;

    std::unique_ptr<ViscousStep> viscousStep;

    ViscousStepTest() {
        // Initialize 10x10x10 grid with spacing 1.0
        grid = std::make_shared<const Grid>(10, 10, 10, 1.0, 1.0, 1.0, *g_mpi);

        // Setup SimulationData context
        data_.grid = grid;
        data_.dt = dt;
        data_.nu = nu;
        data_.currTime = 0.0;
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

        viscousStep = std::make_unique<ViscousStep>(*g_mpi, data_);
        viscousStep->setup();
    }

    void callComputeXi() {
        viscousStep->computeXi();
    }

    double getXi(Axis axis, long i, long j, long k) {
        return viscousStep->xi(axis, i, j, k);
    }
};

TEST_F(ViscousStepTest, ConstructorAndInit) {
    ASSERT_NE(viscousStep, nullptr);
}

TEST_F(ViscousStepTest, RunDoesNotCrash) {
    // Triggers LinearSys::fillSystemVelocity which calls analytical BCs
    ASSERT_NO_THROW(viscousStep->run());
}

TEST_F(ViscousStepTest, ComputeXi_Correctness) {
    // We set up a scenario where the inlined 'g' calculation and 'xi' calculation are predictable.
    // Formula: xi = u + (dt / beta) * g
    // where g = f - gradP - (nu/k)*u + nu * Laplacian
    // and beta = 1 + (dt * nu / 2) * inv_k

    // 1. Force F = 10.0
    auto funcF = [](double, double, double, double) { return 10.0; };
    data_.f.setup(grid, funcF, funcF, funcF);
    data_.f.populate(0.0);

    // 2. Predictor P = 0 => gradP = 0
    auto funcZero = [](double, double, double, double) { return 0.0; };
    data_.predictor.setup(grid, funcZero);
    data_.predictor.populate(0.0);

    // 3. Velocity U = 1.0 (constant)
    auto funcOne = [](double, double, double, double) { return 1.0; };
    data_.u.setup(grid, funcOne, funcOne, funcOne);
    data_.u.populate(0.0);

    // 4. Eta, Zeta = 0 (and U is constant) => Laplacian terms are 0
    data_.eta.setup(grid, funcZero, funcZero, funcZero);
    data_.eta.populate(0.0);
    data_.zeta.setup(grid, funcZero, funcZero, funcZero);
    data_.zeta.populate(0.0);

    // 5. Permeability inv_k = 1.0
    data_.inv_k.setup(grid, funcOne); // Explicitly ensure it is 1.0
    data_.inv_k.populate(0.0);

    // 6. [IMPORTANT] Update BCs to match U=1.0 to avoid ghost cell artifacts
    data_.bcu = funcOne;
    data_.bcv = funcOne;
    data_.bcw = funcOne;

    // --- Execute ---
    callComputeXi();

    // --- Verification ---
    const long i = (long) grid->Nx / 2;
    const long j = (long) grid->Ny / 2;
    const long k_ = (long) grid->Nz / 2;

    // Manual Calculation:
    double u_val = 1.0;
    double f_val = 10.0;
    double gradP_val = 0.0;
    double inv_k = 1.0;
    double laplacian = 0.0;

    // g = 10.0 - 0.0 - (1.0 * 1.0 * 1.0) + 0.0 = 9.0
    double g_val = f_val - gradP_val - nu * u_val * inv_k + nu * laplacian;

    // beta = 1.0 + (0.1 * 1.0 * 0.5 * 1.0) = 1.05
    double beta = 1.0 + (dt * nu * 0.5 * inv_k);

    // xi = 1.0 + (0.1 / 1.05) * 9.0
    double expected_xi = u_val + (dt / beta) * g_val;

    double xi_x_computed = getXi(Axis::X, i, j, k_);
    EXPECT_NEAR(xi_x_computed, expected_xi, 1e-9);
}

class ViscousStepRobustnessTest : public ::testing::Test {
protected:
    const size_t N = 10;
    std::shared_ptr<const Grid> grid;
    SimulationData data_;
    std::unique_ptr<ViscousStep> viscousStep;

    ViscousStepRobustnessTest() {
        grid = std::make_shared<const Grid>(N, N, N, 0.1, 0.1, 0.1, *g_mpi);
        data_.grid = grid;
        data_.dt = 0.1;
        data_.nu = 1.0e-3;
        data_.currTime = 0.0;
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

        viscousStep = std::make_unique<ViscousStep>(*g_mpi, data_);
        viscousStep->setup();
    }

    void checkFieldFinite(const Field &field, const std::string &fieldName) {
        for (long k = 0; k < grid->Nz; ++k) {
            for (long j = 0; j < grid->Ny; ++j) {
                for (long i = 0; i < grid->Nx; ++i) {
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