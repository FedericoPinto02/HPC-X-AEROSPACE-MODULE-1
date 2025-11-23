#include "core/Fields.hpp"
#include "core/Functions.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <memory>

// ============================================================================
// === FIELD CLASS TESTS
// ============================================================================

class FieldTestFixture : public ::testing::Test {
protected:
    const size_t Nx = 4, Ny = 4, Nz = 4;
    const double dx = 1.0, dy = 1.0, dz = 1.0;

    // Grid is held via shared_ptr to match Field's internal storage requirement
    std::shared_ptr<Grid> grid;
    Field testField;

    void SetUp() override {
        grid = std::make_shared<Grid>(Nx, Ny, Nz, dx, dy, dz);

        // Initialize with zero function by default
        testField.setup(grid, Functions::ZERO);
        testField.populate(0.0);
    }
};

TEST_F(FieldTestFixture, GetGrid_ReturnsCorrectAddress) {
    // Verify that the field holds a reference to the exact same grid object
    EXPECT_EQ(&testField.getGrid(), grid.get());
}

TEST_F(FieldTestFixture, Setup_WithFunctions_ZeroDefault) {
    testField.setup(grid, Functions::ZERO);
    testField.populate(0.0);

    // Check corners
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(testField(Nx - 1, Ny - 1, Nz - 1), 0.0);
}

TEST_F(FieldTestFixture, Setup_WithFunctions_CustomLogic) {
    // Define logic: f(x,y,z) = x + 2y + 3z
    auto myFunc = [](double x, double y, double z, double t) {
        return x + 2.0 * y + 3.0 * z;
    };

    testField.setup(grid, myFunc);
    testField.populate(0.0);

    // (0,0,0) -> pos (0.0, 0.0, 0.0) -> Val = 0
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 0.0);

    // (1,0,0) -> pos (1.0, 0.0, 0.0) -> Val = 1.0
    EXPECT_DOUBLE_EQ(testField(1, 0, 0), 1.0);

    // (0,1,0) -> pos (0.0, 1.0, 0.0) -> Val = 2.0
    EXPECT_DOUBLE_EQ(testField(0, 1, 0), 2.0);
}

TEST_F(FieldTestFixture, Populate_TimeDependence) {
    // Define time-dependent logic: f(x,y,z,t) = t
    testField.setup(grid, [](double, double, double, double t) { return t; });

    // Populate at t=10.0
    testField.populate(10.0);
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 10.0);

    // Populate at t=25.5
    testField.populate(25.5);
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 25.5);
}

TEST_F(FieldTestFixture, Access_StandardIndexing) {
    // Map grid coordinates directly to value to verify indexing
    testField.setup(grid, [](double x, double, double, double) { return x; });
    testField.populate(0.0);

    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(testField(1, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(testField(Nx - 1, 0, 0), static_cast<double>(Nx - 1));
}

TEST_F(FieldTestFixture, Access_ValueWithOffset) {
    // f(x, y, z) = x + 100y + 10000z
    auto distFunc = [](double x, double y, double z, double) {
        return x + 100.0 * y + 10000.0 * z;
    };

    testField.setup(grid, distFunc);
    testField.populate(0.0);

    double centerVal = testField(2, 2, 2);

    // Check neighbor in X+1 direction (should increase value by 1.0)
    EXPECT_DOUBLE_EQ(testField.valueWithOffset(2, 2, 2, Axis::X, 1), centerVal + 1.0);

    // Check neighbor in Y-1 direction (should decrease value by 100.0)
    EXPECT_DOUBLE_EQ(testField.valueWithOffset(2, 2, 2, Axis::Y, -1), centerVal - 100.0);
}

TEST_F(FieldTestFixture, Add_Scalar) {
    testField.reset(10.0);
    testField.add(5.0);
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 15.0);
}

TEST_F(FieldTestFixture, Add_Field_MismatchThrows) {
    // Create a grid with different dimensions
    auto wrongGrid = std::make_shared<Grid>(Nx + 1, Ny, Nz, 1.0, 1.0, 1.0);
    Field wrongField;

    wrongField.setup(wrongGrid, Functions::ZERO);
    wrongField.populate(0.0);

    // Expect error when adding fields from different grids
    EXPECT_THROW(testField.add(wrongField), std::invalid_argument);
}

// ============================================================================
// === VECTOR FIELD TESTS
// ============================================================================

class VectorFieldTestFixture : public ::testing::Test {
protected:
    const size_t Nx = 2, Ny = 2, Nz = 2;
    std::shared_ptr<Grid> grid;
    VectorField vField;

    void SetUp() override {
        grid = std::make_shared<Grid>(Nx, Ny, Nz, 1.0, 1.0, 1.0);
        vField.setup(grid);
    }
};

TEST_F(VectorFieldTestFixture, Setup_DistinctFunctionsPerComponent) {
    // Define distinct functions for X, Y, Z components
    auto fX = [](double, double, double, double) { return 1.0; };
    auto fY = [](double, double, double, double) { return 2.0; };
    auto fZ = [](double, double, double, double) { return 3.0; };

    vField.setup(grid, fX, fY, fZ);
    vField.populate(0.0);

    EXPECT_DOUBLE_EQ(vField(Axis::X)(0, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(vField(Axis::Y)(0, 0, 0), 2.0);
    EXPECT_DOUBLE_EQ(vField(Axis::Z)(0, 0, 0), 3.0);
}

TEST_F(VectorFieldTestFixture, Add_VectorField) {
    auto ones = [](double, double, double, double) { return 1.0; };
    vField.setup(grid, ones, ones, ones);
    vField.populate();

    VectorField v2;
    auto twos = [](double, double, double, double) { return 2.0; };
    v2.setup(grid, twos, twos, twos);
    v2.populate();

    // Add v2 to vField: (1,1,1) + (2,2,2) = (3,3,3)
    vField.add(v2);

    EXPECT_DOUBLE_EQ(vField(Axis::X)(0, 0, 0), 3.0);
    EXPECT_DOUBLE_EQ(vField(Axis::Y)(0, 0, 0), 3.0);
    EXPECT_DOUBLE_EQ(vField(Axis::Z)(0, 0, 0), 3.0);
}