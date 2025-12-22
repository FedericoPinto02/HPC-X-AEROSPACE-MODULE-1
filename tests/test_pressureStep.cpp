#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>

#include "MpiEnvFixture.hpp"

#include "simulation/pressureStep.hpp"
#include "simulation/SimulationContext.hpp"
#include "core/Fields.hpp"


class PressureStepTest : public ::testing::Test {
protected:
    const double dt = 0.1;

    // Grid managed via shared_ptr to match Field's internal storage requirement
    std::shared_ptr<const Grid> grid;

    // Context objects
    SimulationData data_;

    std::unique_ptr<PressureStep> pressureStep;

    PressureStepTest() {
        assert(g_mpi != nullptr);

        // Initialize 10x10x10 grid with spacing 1.0
        grid = std::make_shared<const Grid>(10, 10, 10, 1.0, 1.0, 1.0, *g_mpi);

        // Setup SimulationData context
        data_.grid = grid;
        data_.dt = dt;
    }

    void SetUp() override {
        // Define Pressure: p(x,y,z) = 0.0
        data_.p.setup(grid, ZERO_FUNC);
        data_.p.populate(0.0);
        data_.predictor.setup(grid, ZERO_FUNC);
        data_.predictor.populate(0.0);

        // Define Velocity: u(x,y,z) = [sin(x), cos(y), sin(z)]
        auto funcX = [](double x, double, double, double) { return std::sin(x); };
        auto funcY = [](double, double y, double, double) { return std::cos(y); };
        auto funcZ = [](double, double, double z, double) { return std::sin(z); };

        data_.u.setup(grid, funcX, funcY, funcZ);
        data_.u.populate(0.0);
        data_.bcu = funcX;
        data_.bcv = funcY;
        data_.bcw = funcZ;

        // Create the object under test
        pressureStep = std::make_unique<PressureStep>(*g_mpi, data_);
        pressureStep->setup();
    }

    double getDivU(long i, long j, long k) {
        return pressureStep->divU(i, j, k);
    }

    double getPressure(long i, long j, long k) {
        return data_.p(i, j, k);
    }
};

TEST_F(PressureStepTest, ConstructorAndInit) {
    ASSERT_NE(pressureStep, nullptr);
}

TEST_F(PressureStepTest, RunDoesNotCrash) {
    ASSERT_NO_THROW(pressureStep->run());
}

TEST_F(PressureStepTest, Run_CalculatesDivU_Correctness) {
    ASSERT_NO_THROW(pressureStep->run());

    const long i = (long) grid->Nx / 2;
    const long j = (long) grid->Ny / 2;
    const long k = (long) grid->Nz / 2;
    const double dx = grid->dx, dy = grid->dy, dz = grid->dz;

    // 1. X-Term: U is staggered in X.
    // We want (u[i] - u[i-1])/dx.
    // u[i] is at (i+0.5)dx, u[i-1] is at (i-0.5)dx
    double x_curr = grid->to_x(i, GridStaggering::FACE_CENTERED, Axis::X);
    double x_prev = grid->to_x(i - 1, GridStaggering::FACE_CENTERED, Axis::X);
    double du_dx = (std::sin(x_curr) - std::sin(x_prev)) / dx;

    // 2. Y-Term: V is staggered in Y.
    double y_curr = grid->to_y(j, GridStaggering::FACE_CENTERED, Axis::Y);
    double y_prev = grid->to_y(j - 1, GridStaggering::FACE_CENTERED, Axis::Y);
    double du_dy = (std::cos(y_curr) - std::cos(y_prev)) / dy;

    // 3. Z-Term: W is staggered in Z.
    double z_curr = grid->to_z(k, GridStaggering::FACE_CENTERED, Axis::Z);
    double z_prev = grid->to_z(k - 1, GridStaggering::FACE_CENTERED, Axis::Z);
    double du_dz = (std::sin(z_curr) - std::sin(z_prev)) / dz;

    double expected_divU = du_dx + du_dy + du_dz;

    double divU_computed = getDivU(i, j, k);
    EXPECT_NEAR(divU_computed, expected_divU, 1e-9);
}

TEST_F(PressureStepTest, Run_WithZeroDivergence_PressureIsUnchanged) {
    const long i = (long) grid->Nx / 2;
    const long j = (long) grid->Ny / 2;
    const long k = (long) grid->Nz / 2;

    // Setup U = 0, P = 1.0
    data_.u.setup(grid, ZERO_FUNC, ZERO_FUNC, ZERO_FUNC);
    data_.u.populate(0.0);

    data_.bcu = ZERO_FUNC;
    data_.bcv = ZERO_FUNC;
    data_.bcw = ZERO_FUNC;

    data_.p.setup(grid, [](double, double, double, double) { return 1.0; });
    data_.p.populate(0.0);

    // Reconstruct with new data
    pressureStep = std::make_unique<PressureStep>(*g_mpi, data_);
    pressureStep->setup();

    pressureStep->run();

    // If div(u) = 0, pcr should be 0, so p_final == p_initial == 1.0
    double p_final = getPressure(i, j, k);
    EXPECT_NEAR(p_final, 1.0, 1e-9);

    double divU_computed = getDivU(i, j, k);
    EXPECT_NEAR(divU_computed, 0.0, 1e-9);
}

class PressureStepRobustnessTest : public ::testing::Test {
protected:
    const size_t N = 10;
    std::shared_ptr<Grid> grid;
    SimulationData data_;
    std::unique_ptr<PressureStep> pressureStep;

    PressureStepRobustnessTest() {
        assert(g_mpi != nullptr);

        grid = std::make_shared<Grid>(N, N, N, 0.1, 0.1, 0.1, *g_mpi);
        data_.grid = grid;
        data_.dt = 0.1;
    }

    void SetUp() override {
        // Define Pressure: p = x*y + z
        auto funcP = [](double x, double y, double z, double) {
            return x * y + z;
        };
        data_.p.setup(grid, funcP);
        data_.p.populate(0.0);

        data_.predictor.setup(grid, ZERO_FUNC);
        data_.predictor.populate(0.0);

        // Define Velocity: u = [x^2, y^2, z^2]
        auto funcU_x = [](double x, double, double, double) { return x * x; };
        auto funcU_y = [](double, double y, double, double) { return y * y; };
        auto funcU_z = [](double, double, double z, double) { return z * z; };

        data_.u.setup(grid, funcU_x, funcU_y, funcU_z);
        data_.u.populate(0.0);
        data_.bcu = funcU_x;
        data_.bcv = funcU_y;
        data_.bcw = funcU_z;

        pressureStep = std::make_unique<PressureStep>(*g_mpi, data_);
        pressureStep->setup();
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

TEST_F(PressureStepRobustnessTest, Run_WithPolynomialFields_ProducesFiniteOutput) {
    ASSERT_NO_THROW(pressureStep->run());
    checkFieldFinite(data_.p, "data_.p");
}
