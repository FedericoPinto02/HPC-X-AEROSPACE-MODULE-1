#include "core/Fields.hpp"
#include "core/Functions.hpp"
#include <cmath>
#include <gtest/gtest.h>

class FieldTestFixture : public ::testing::Test {
protected:
    const size_t Nx = 4, Ny = 4, Nz = 4;
    const double dx = 1.0, dy = 1.0, dz = 1.0;

    const Grid grid;
    Field testField;

    FieldTestFixture() : grid(Nx, Ny, Nz, dx, dy, dz) {}

    void SetUp() override {
        testField.setup(grid, Functions::ZERO);
        testField.populate(0.0);
    }
};

TEST_F(FieldTestFixture, GetGrid_ReturnsCorrectAddress) {
    EXPECT_EQ(&testField.getGrid(), &grid);
}

TEST_F(FieldTestFixture, Setup_WithFunctions_ZeroDefault) {
    testField.setup(grid, Functions::ZERO);
    testField.populate(0.0);

    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(testField(Nx - 1, Ny - 1, Nz - 1), 0.0);
}

TEST_F(FieldTestFixture, Setup_WithFunctions_CustomLogic) {
    // f(x,y,z) = x + 2y + 3z
    // Grid dx=dy=dz=1.0. Cell Centered offset is 0.0.
    auto myFunc = [](double x, double y, double z, double t) {
        return x + 2.0 * y + 3.0 * z;
    };

    testField.setup(grid, myFunc);
    testField.populate(0.0);

    // (0,0,0) -> pos (0.0, 0.0, 0.0) -> Val = 0 + 0 + 0 = 0.0
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 0.0);

    // (1,0,0) -> pos (1.0, 0.0, 0.0) -> Val = 1.0 + 0 + 0 = 1.0
    EXPECT_DOUBLE_EQ(testField(1, 0, 0), 1.0);

    // (0,1,0) -> pos (0.0, 1.0, 0.0) -> Val = 0 + 2.0 + 0 = 2.0
    EXPECT_DOUBLE_EQ(testField(0, 1, 0), 2.0);
}

TEST_F(FieldTestFixture, Populate_TimeDependence) {
    // f(x,y,z,t) = t
    testField.setup(grid, [](double, double, double, double t) { return t; });

    testField.populate(10.0);
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 10.0);

    testField.populate(25.5);
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 25.5);
}

TEST_F(FieldTestFixture, Access_StandardIndexing) {
    // Map grid coordinates directly to value to verify indexing
    // With Cell Centered offset 0.0: x = i * dx
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

    // Center point indices (2, 2, 2)
    double centerVal = testField(2, 2, 2);

    // Neighbor in X+1 -> physical x increases by 1.0 -> value increases by 1.0
    EXPECT_DOUBLE_EQ(testField.valueWithOffset(2, 2, 2, Axis::X, 1), centerVal + 1.0);

    // Neighbor in Y-1 -> physical y decreases by 1.0 -> value decreases by 100.0
    EXPECT_DOUBLE_EQ(testField.valueWithOffset(2, 2, 2, Axis::Y, -1), centerVal - 100.0);
}

TEST_F(FieldTestFixture, Add_Scalar) {
    testField.reset(10.0);
    testField.add(5.0);
    EXPECT_DOUBLE_EQ(testField(0, 0, 0), 15.0);
}

TEST_F(FieldTestFixture, Add_Field_MismatchThrows) {
    Grid wrongGrid(Nx + 1, Ny, Nz, 1.0, 1.0, 1.0);
    Field wrongField;

    wrongField.setup(wrongGrid, Functions::ZERO);
    wrongField.populate(0.0);

    EXPECT_THROW(testField.add(wrongField), std::invalid_argument);
}

class VectorFieldTestFixture : public ::testing::Test {
protected:
    const size_t Nx = 2, Ny = 2, Nz = 2;
    const Grid grid;
    VectorField vField;

    VectorFieldTestFixture() : grid(Nx, Ny, Nz, 1.0, 1.0, 1.0) {}

    void SetUp() override {
        vField.setup(grid);
    }
};

TEST_F(VectorFieldTestFixture, Setup_DistinctFunctionsPerComponent) {
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

    vField.add(v2);

    EXPECT_DOUBLE_EQ(vField(Axis::X)(0, 0, 0), 3.0);
    EXPECT_DOUBLE_EQ(vField(Axis::Y)(0, 0, 0), 3.0);
    EXPECT_DOUBLE_EQ(vField(Axis::Z)(0, 0, 0), 3.0);
}