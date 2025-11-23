#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>
#include <functional>

#include "simulation/viscousStep.hpp"
#include "simulation/SimulationContext.hpp"
#include "core/Fields.hpp"
#include "core/Functions.hpp"
#include "numerics/derivatives.hpp"

class ViscousStepTest : public ::testing::Test {
protected:
    const double dt = 0.1;
    const double nu = 1.0;

    // Stack-allocated Grid
    Grid grid;

    // Context
    SimulationData data_;
    ParallelizationSettings parallel_;

    std::unique_ptr<ViscousStep> viscousStep;

    ViscousStepTest()
            : grid(10, 10, 10, 1.0, 1.0, 1.0),
              data_{ grid }
    {
        data_.dt = dt;
        data_.nu = nu;
        data_.currTime = 0.0; // Ensure time is initialized
    }

    void SetUp() override {
        // k = 1.0 everywhere
        data_.k.setup(grid, [](double,double,double,double){ return 1.0; });
        data_.k.populate(0.0);

        // f = [0,0,0]
        data_.f.setup(grid, Functions::ZERO, Functions::ZERO, Functions::ZERO);
        data_.f.populate(0.0);

        // 1. Pressure: p(x,y,z) = x + 2y + 3z
        auto funcP = [](double x, double y, double z, double) {
            return x + 2.0 * y + 3.0 * z;
        };
        data_.p.setup(grid, funcP);
        data_.p.populate(0.0);

        // 2. Velocity: u = [sin(x), sin(y), sin(z)]
        auto sinX = [](double x, double, double, double) { return std::sin(x); };
        auto sinY = [](double, double y, double, double) { return std::sin(y); };
        auto sinZ = [](double, double, double z, double) { return std::sin(z); };
        data_.u.setup(grid, sinX, sinY, sinZ);
        data_.u.populate(0.0);

        // 3. Eta: eta = [cos(x), cos(y), cos(z)]
        auto cosX = [](double x, double, double, double) { return std::cos(x); };
        auto cosY = [](double, double y, double, double) { return std::cos(y); };
        auto cosZ = [](double, double, double z, double) { return std::cos(z); };
        data_.eta.setup(grid, cosX, cosY, cosZ);
        data_.eta.populate(0.0);

        // 4. Zeta: zeta = [x^2, y^2, z^2]
        auto sqX = [](double x, double, double, double) { return x * x; };
        auto sqY = [](double, double y, double, double) { return y * y; };
        auto sqZ = [](double, double, double z, double) { return z * z; };
        data_.zeta.setup(grid, sqX, sqY, sqZ);
        data_.zeta.populate(0.0);

        // Boundary Fields (Discrete)
        data_.uBoundNew.setup(grid, sinX, sinY, sinZ);
        data_.uBoundNew.populate(0.0);
        data_.uBoundOld.setup(grid, sinX, sinY, sinZ);
        data_.uBoundOld.populate(0.0);

        // --- FIX: Initialize Analytical Boundary Functions ---
        // These are required by LinearSys::fillSystemVelocity (lines 100+)
        // matching the velocity definitions above.
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
    // This call triggers LinearSys::fillSystemVelocity, which calls data_.bcu/v/w
    ASSERT_NO_THROW(viscousStep->run());
}

TEST_F(ViscousStepTest, ComputeG_Correctness) {
    const size_t i = 5, j = 5, k_ = 5;
    const double x = i * grid.dx;

    ASSERT_NO_THROW(callComputeG());

    double f_term = 0.0;
    double gradP_term = 1.0;
    double u_term = (nu / 1.0) * std::sin(x);

    double dxx_term = std::cos(x + 1.0) - 2.0 * std::cos(x) + std::cos(x - 1.0);
    double dyy_term = 0.0;
    double dzz_term = 0.0;

    double deriv_term = nu * (dxx_term + dyy_term + dzz_term);

    double expected_g_x = f_term - gradP_term - u_term * 0.5 + deriv_term * 0.5;

    double g_x_computed = getG(Axis::X, i, j, k_);
    EXPECT_NEAR(g_x_computed, expected_g_x, 1e-9);
}

TEST_F(ViscousStepTest, ComputeXi_Correctness) {
    const size_t i = 5, j = 5, k_ = 5;
    const double x = i * grid.dx;
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
    Grid grid;
    SimulationData data_;
    ParallelizationSettings parallel_;
    std::unique_ptr<ViscousStep> viscousStep;

    ViscousStepRobustnessTest()
            : grid(N, N, N, 0.1, 0.1, 0.1),
              data_{ grid }
    {
        data_.dt = 0.1;
        data_.nu = 1.0e-3;
        data_.currTime = 0.0;
    }

    void SetUp() override {
        // k = 1.0
        data_.k.setup(grid, [](double,double,double,double){ return 1.0; });
        data_.k.populate(0.0);

        // f = 10.0
        auto funcF = [](double,double,double,double){ return 10.0; };
        data_.f.setup(grid, funcF, funcF, funcF);
        data_.f.populate(0.0);

        // p = x*y + z
        data_.p.setup(grid, [](double x, double y, double z, double){ return x * y + z; });
        data_.p.populate(0.0);

        // u = [x^2, y^2, z^2]
        auto sqX = [](double x, double, double, double){ return x*x; };
        auto sqY = [](double, double y, double, double){ return y*y; };
        auto sqZ = [](double, double, double z, double){ return z*z; };
        data_.u.setup(grid, sqX, sqY, sqZ);
        data_.u.populate(0.0);

        // eta = [y^2, z^2, x^2]
        data_.eta.setup(grid, sqY, sqZ, sqX);
        data_.eta.populate(0.0);

        // zeta = [z^2, x^2, y^2]
        data_.zeta.setup(grid, sqZ, sqX, sqY);
        data_.zeta.populate(0.0);

        // BCs
        data_.uBoundNew.setup(grid, sqX, sqY, sqZ);
        data_.uBoundNew.populate(0.0);
        data_.uBoundOld.setup(grid, sqX, sqY, sqZ);
        data_.uBoundOld.populate(0.0);

        // --- FIX: Initialize Analytical Boundary Functions ---
        data_.bcu = sqX;
        data_.bcv = sqY;
        data_.bcw = sqZ;

        viscousStep = std::make_unique<ViscousStep>(data_, parallel_);
    }

    void checkFieldFinite(const Field& field, const std::string& fieldName) {
        for (size_t k = 0; k < grid.Nz; ++k) {
            for (size_t j = 0; j < grid.Ny; ++j) {
                for (size_t i = 0; i < grid.Nx; ++i) {
                    ASSERT_TRUE(std::isfinite(field(i, j, k)))
                                                << "Found NaN or Inf in " << fieldName
                                                << " at index (" << i << "," << j << "," << k << ")";
                }
            }
        }
    }
};

TEST_F(ViscousStepRobustnessTest, Run_WithPolynomialFields_ProducesFiniteOutput)
{
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