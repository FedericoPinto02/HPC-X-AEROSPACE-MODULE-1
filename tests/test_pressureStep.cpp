#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>
#include <functional>

#include "simulation/pressureStep.hpp"
#include "simulation/SimulationContext.hpp"
#include "core/Fields.hpp"
#include "core/Functions.hpp"
#include "numerics/derivatives.hpp"

class PressureStepTest : public ::testing::Test {
protected:
    const double dt = 0.1;

    // Grid managed via shared_ptr to match Field's internal storage requirement
    std::shared_ptr<Grid> grid;

    // Context objects
    SimulationData data_;
    ParallelizationSettings parallel_;

    std::unique_ptr<PressureStep> pressureStep;

    PressureStepTest() {
        // Initialize 10x10x10 grid with spacing 1.0
        grid = std::make_shared<Grid>(10, 10, 10, 1.0, 1.0, 1.0);

        // Setup SimulationData context
        data_.grid = grid;
        data_.dt = dt;

        parallel_.schurDomains = 1;
    }

    void SetUp() override {
        // Define Pressure: p(x,y,z) = 0.0
        data_.p.setup(grid, Functions::ZERO);
        data_.p.populate(0.0);

        // Define Velocity: u(x,y,z) = [sin(x), cos(y), sin(z)]
        auto funcX = [](double x, double, double, double) { return std::sin(x); };
        auto funcY = [](double, double y, double, double) { return std::cos(y); };
        auto funcZ = [](double, double, double z, double) { return std::sin(z); };

        data_.u.setup(grid, funcX, funcY, funcZ);
        data_.u.populate(0.0);

        // Create the object under test
        pressureStep = std::make_unique<PressureStep>(data_, parallel_);
    }

    double getDivU(size_t i, size_t j, size_t k) {
        return pressureStep->divU(i, j, k);
    }

    double getPressure(size_t i, size_t j, size_t k) {
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

    const size_t i = 5, j = 5, k = 5;
    const double x = i * grid->dx, y = j * grid->dy, z = k * grid->dz;
    const double dx = grid->dx, dy = grid->dy, dz = grid->dz;

    // The Derivatives class uses Backward Difference for divergence
    // u.x = sin(x) -> d/dx ~ (sin(x) - sin(x-dx)) / dx
    double du_dx = (std::sin(x) - std::sin(x - dx)) / dx;

    // u.y = cos(y) -> d/dy ~ (cos(y) - cos(y-dy)) / dy
    double du_dy = (std::cos(y) - std::cos(y - dy)) / dy;

    // u.z = sin(z) -> d/dz ~ (sin(z) - sin(z-dz)) / dz
    double du_dz = (std::sin(z) - std::sin(z - dz)) / dz;

    double expected_divU = du_dx + du_dy + du_dz;

    double divU_computed = getDivU(i, j, k);
    EXPECT_NEAR(divU_computed, expected_divU, 1e-9);
}

TEST_F(PressureStepTest, Run_WithZeroDivergence_PressureIsUnchanged) {
    // Setup U = 0, P = 1.0
    data_.u.setup(grid, Functions::ZERO, Functions::ZERO, Functions::ZERO);
    data_.u.populate(0.0);

    data_.p.setup(grid, [](double, double, double, double) { return 1.0; });
    data_.p.populate(0.0);

    // Reconstruct with new data
    pressureStep = std::make_unique<PressureStep>(data_, parallel_);

    pressureStep->run();

    // If div(u) = 0, pcr should be 0, so p_final == p_initial == 1.0
    double p_final = getPressure(5, 5, 5);
    EXPECT_NEAR(p_final, 1.0, 1e-9);

    double divU_computed = getDivU(5, 5, 5);
    EXPECT_NEAR(divU_computed, 0.0, 1e-9);
}

class PressureStepRobustnessTest : public ::testing::Test {
protected:
    const size_t N = 10;
    std::shared_ptr<Grid> grid;
    SimulationData data_;
    ParallelizationSettings parallel_;
    std::unique_ptr<PressureStep> pressureStep;

    PressureStepRobustnessTest() {
        grid = std::make_shared<Grid>(N, N, N, 0.1, 0.1, 0.1);
        data_.grid = grid;
        data_.dt = 0.1;
        parallel_.schurDomains = 1;
    }

    void SetUp() override {
        // Define Pressure: p = x*y + z
        auto funcP = [](double x, double y, double z, double) {
            return x * y + z;
        };
        data_.p.setup(grid, funcP);
        data_.p.populate(0.0);

        // Define Velocity: u = [x^2, y^2, z^2]
        auto funcU_x = [](double x, double, double, double) { return x * x; };
        auto funcU_y = [](double, double y, double, double) { return y * y; };
        auto funcU_z = [](double, double, double z, double) { return z * z; };

        data_.u.setup(grid, funcU_x, funcU_y, funcU_z);
        data_.u.populate(0.0);

        pressureStep = std::make_unique<PressureStep>(data_, parallel_);
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

TEST_F(PressureStepRobustnessTest, Run_WithPolynomialFields_ProducesFiniteOutput) {
    ASSERT_NO_THROW(pressureStep->run());
    checkFieldFinite(data_.p, "data_.p");
}