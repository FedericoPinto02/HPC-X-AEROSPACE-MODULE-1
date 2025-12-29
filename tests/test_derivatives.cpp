#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <functional>
#include <memory>

#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"
#include "core/Grid.hpp"

// dummy BC function (returns 0.0)
const Func dummyBC = 
[](double x, double y, double z, double t) { return 0.0; };
double simTime = 0.0;

class DerivativesTest : public ::testing::Test {
protected:
    // Constants for grid dimensions and spacing
    const size_t Nx = 10, Ny = 10, Nz = 10;
    const double dx = 0.5, dy = 0.5, dz = 0.5;

    // Grid managed via shared_ptr to match Field architecture
    std::shared_ptr<Grid> grid;
    Derivatives deriv;

    void SetUp() override {
        grid = std::make_shared<Grid>(Nx, Ny, Nz, dx, dy, dz);
    }

    // Helper to quickly setup and populate a field with a mathematical function
    Field createField(Func func) {
        Field f;
        f.setup(grid, func);
        f.populate(0.0);
        return f;
    }

    // Helper to setup an empty field
    Field createEmptyField() {
        Field f;
        f.setup(grid);
        return f;
    }

    // Helper to setup an empty VectorField
    VectorField createEmptyVectorField() {
        VectorField vf;
        vf.setup(grid);
        return vf;
    }
};

// ============================================================================
// === SCALAR FIELD: First-Order Forward Differences (Dx, Dy, Dz)
// ============================================================================

TEST_F(DerivativesTest, Scalar_Dx_Forward_Accuracy) {
    // f(x) = x^2
    // Forward Difference: (f(x+h) - f(x))/h = (x^2 + 2xh + h^2 - x^2)/h = 2x + h
    Field field = createField([](double x, double, double, double) { return x * x; });
    Field dxField = createEmptyField();

    deriv.computeDx_fwd(field, dxField);

    const double tol = 1e-6;
    for (size_t i = 0; i < Nx - 1; ++i) { // Valid up to N-2
        double x = i * dx;
        double expected = 2.0 * x + dx;
        EXPECT_NEAR(dxField(i, 0, 0), expected, tol);
    }
    // Boundary check (implementation dependent, usually 0)
    EXPECT_DOUBLE_EQ(dxField(Nx - 1, 0, 0), 0.0);
}

TEST_F(DerivativesTest, Scalar_Dy_Forward_Accuracy) {
    // f(y) = y^2
    // Expected: 2y + dy
    Field field = createField([](double, double y, double, double) { return y * y; });
    Field dyField = createEmptyField();

    deriv.computeDy_fwd(field, dyField);

    const double tol = 1e-6;
    for (size_t j = 0; j < Ny - 1; ++j) {
        double y = j * dy;
        double expected = 2.0 * y + dy;
        EXPECT_NEAR(dyField(0, j, 0), expected, tol);
    }
    EXPECT_DOUBLE_EQ(dyField(0, Ny - 1, 0), 0.0);
}

TEST_F(DerivativesTest, Scalar_Dz_Forward_Accuracy) {
    // f(z) = z^2
    // Expected: 2z + dz
    Field field = createField([](double, double, double z, double) { return z * z; });
    Field dzField = createEmptyField();

    deriv.computeDz_fwd(field, dzField);

    const double tol = 1e-6;
    for (size_t k = 0; k < Nz - 1; ++k) {
        double z = k * dz;
        double expected = 2.0 * z + dz;
        EXPECT_NEAR(dzField(0, 0, k), expected, tol);
    }
    EXPECT_DOUBLE_EQ(dzField(0, 0, Nz - 1), 0.0);
}

// ============================================================================
// === SCALAR FIELD: First-Order Backward Differences (Dx, Dy, Dz)
// ============================================================================

TEST_F(DerivativesTest, Scalar_Dx_Backward_Accuracy) {
    // f(x) = x^2
    // Backward Difference: (f(x) - f(x-h))/h = (x^2 - (x^2 - 2xh + h^2))/h = 2x - h
    Field field = createField([](double x, double, double, double) { return x * x; });
    Field dxField = createEmptyField();

    // Define Boundary Condition: Value of function at boundary x=0 is 0.0
    Func bcFunc = [](double, double, double, double) { return 0.0; };
    double time = 0.0;

    deriv.computeDx_bwd(field, dxField, bcFunc, time);

    const double tol = 1e-6;
    for (size_t i = 1; i < Nx; ++i) { // Valid starting from 1
        double x = i * dx;
        double expected = 2.0 * x - dx;
        EXPECT_NEAR(dxField(i, 0, 0), expected, tol);
    }
}

TEST_F(DerivativesTest, Scalar_Dy_Backward_Accuracy) {
    // f(y) = y^2
    // Expected: 2y - dy
    Field field = createField([](double, double y, double, double) { return y * y; });
    Field dyField = createEmptyField();

    // Define Boundary Condition: Value of function at boundary y=0 is 0.0
    Func bcFunc = [](double, double, double, double) { return 0.0; };
    double time = 0.0;

    deriv.computeDy_bwd(field, dyField, bcFunc, time);

    const double tol = 1e-6;
    for (size_t j = 1; j < Ny; ++j) {
        double y = j * dy;
        double expected = 2.0 * y - dy;
        EXPECT_NEAR(dyField(0, j, 0), expected, tol);
    }
}

TEST_F(DerivativesTest, Scalar_Dz_Backward_Accuracy) {
    // f(z) = z^2
    // Expected: 2z - dz
    Field field = createField([](double, double, double z, double) { return z * z; });
    Field dzField = createEmptyField();

    // Define Boundary Condition: Value of function at boundary z=0 is 0.0
    Func bcFunc = [](double, double, double, double) { return 0.0; };
    double time = 0.0;

    deriv.computeDz_bwd(field, dzField, bcFunc, time);

    const double tol = 1e-6;
    for (size_t k = 1; k < Nz; ++k) {
        double z = k * dz;
        double expected = 2.0 * z - dz;
        EXPECT_NEAR(dzField(0, 0, k), expected, tol);
    }
}

// ============================================================================
// === SCALAR FIELD: Second-Order Central Differences (Dxx, Dyy, Dzz)
// ============================================================================

TEST_F(DerivativesTest, Scalar_Dxx_Accuracy) {
    // f(x) = x^3. Exact Dxx = 6x.
    // Central difference is exact for cubic polynomials.
    Field field = createField([](double x, double, double, double) { return x * x * x; });
    Field dxx = createEmptyField();

    deriv.computeDxx(field, dxx, dummyBC, simTime, Axis::X);

    const double tol = 1e-8;
    for (size_t i = 1; i < Nx - 1; ++i) { // Valid interior [1, N-2]
        double x = i * dx;
        EXPECT_NEAR(dxx(i, 0, 0), 6.0 * x, tol);
    }
}

TEST_F(DerivativesTest, Scalar_Dyy_Accuracy) {
    // f(y) = y^3. Exact Dyy = 6y.
    Field field = createField([](double, double y, double, double) { return y * y * y; });
    Field dyy = createEmptyField();

    deriv.computeDyy(field, dyy, dummyBC, simTime, Axis::Y);

    const double tol = 1e-8;
    for (size_t j = 1; j < Ny - 1; ++j) {
        double y = j * dy;
        EXPECT_NEAR(dyy(0, j, 0), 6.0 * y, tol);
    }
}

TEST_F(DerivativesTest, Scalar_Dzz_Accuracy) {
    // f(z) = z^3. Exact Dzz = 6z.
    Field field = createField([](double, double, double z, double) { return z * z * z; });
    Field dzz = createEmptyField();

    deriv.computeDzz(field, dzz, dummyBC, simTime, Axis::Z);

    const double tol = 1e-8;
    for (size_t k = 1; k < Nz - 1; ++k) {
        double z = k * dz;
        EXPECT_NEAR(dzz(0, 0, k), 6.0 * z, tol);
    }
}

// ============================================================================
// === VECTOR FIELD: Second-Order Central Differences
// ============================================================================

TEST_F(DerivativesTest, Vector_Dxx_Accuracy) {
    // Verify Dxx operates on ALL components of the vector field.
    // u_x = x^3, u_y = x^3, u_z = x^3
    // Expect dxx = 6x for all components.
    auto func = [](double x, double, double, double) { return x * x * x; };

    VectorField vec;
    vec.setup(grid, func, func, func);
    vec.populate(0.0);

    VectorField result;
    result.setup(grid);

    deriv.computeDxx(vec, result, dummyBC, dummyBC, dummyBC, simTime);

    const double tol = 1e-8;
    for (size_t i = 1; i < Nx - 1; ++i) {
        double x_center = i * dx;
        double x_staggered = (i + 0.5) * dx; // Shift for FACE_CENTERED X-component

        // X component IS staggered in X -> Expected 6 * (x + 0.5dx)
        EXPECT_NEAR(result(Axis::X, i, 0, 0), 6.0 * x_staggered, tol);

        // Y and Z components are NOT staggered in X -> Expected 6 * x
        EXPECT_NEAR(result(Axis::Y, i, 0, 0), 6.0 * x_center, tol);
        EXPECT_NEAR(result(Axis::Z, i, 0, 0), 6.0 * x_center, tol);
    }
}

TEST_F(DerivativesTest, Vector_Dyy_Accuracy) {
    // u = [y^3, y^3, y^3]
    auto func = [](double, double y, double, double) { return y * y * y; };
    VectorField vec;
    vec.setup(grid, func, func, func);
    vec.populate(0.0);

    VectorField result;
    result.setup(grid);

    deriv.computeDyy(vec, result, dummyBC, dummyBC, dummyBC, simTime);

    const double tol = 1e-8;
    for (size_t j = 1; j < Ny - 1; ++j) {
        double y_center = j * dy;
        double y_staggered = (j + 0.5) * dy; // Shift for FACE_CENTERED Y-component

        // X and Z components are NOT staggered in Y -> Expected 6 * y
        EXPECT_NEAR(result(Axis::X, 0, j, 0), 6.0 * y_center, tol);
        EXPECT_NEAR(result(Axis::Z, 0, j, 0), 6.0 * y_center, tol);

        // Y component IS staggered in Y -> Expected 6 * (y + 0.5dy)
        EXPECT_NEAR(result(Axis::Y, 0, j, 0), 6.0 * y_staggered, tol);
    }
}

TEST_F(DerivativesTest, Vector_Dzz_Accuracy) {
    // u = [z^3, z^3, z^3]
    auto func = [](double, double, double z, double) { return z * z * z; };
    VectorField vec;
    vec.setup(grid, func, func, func);
    vec.populate(0.0);

    VectorField result;
    result.setup(grid);

    deriv.computeDzz(vec, result, dummyBC, dummyBC, dummyBC, simTime);

    const double tol = 1e-8;
    for (size_t k = 1; k < Nz - 1; ++k) {
        double z_center = k * dz;
        double z_staggered = (k + 0.5) * dz; // Shift for FACE_CENTERED

        // X and Y components are NOT staggered in Z -> Use Center Z
        EXPECT_NEAR(result(Axis::X, 0, 0, k), 6.0 * z_center, tol);
        EXPECT_NEAR(result(Axis::Y, 0, 0, k), 6.0 * z_center, tol);

        // Z component IS staggered in Z -> Use Staggered Z
        EXPECT_NEAR(result(Axis::Z, 0, 0, k), 6.0 * z_staggered, tol);
    }
}

// ============================================================================
// === COMPOSITE METHODS: Gradient & Divergence
// ============================================================================

TEST_F(DerivativesTest, Gradient_Composite_Test) {
    // f(x,y,z) = x + 2y + 3z
    // Gradient should be [1, 2, 3] everywhere
    Field field = createField([](double x, double y, double z, double) {
        return x + 2.0 * y + 3.0 * z;
    });

    VectorField grad = createEmptyVectorField();
    deriv.computeGradient(field, grad);

    // Check internal points (Forward Difference for linear is exact if grid is uniform)
    const double tol = 1e-8;
    // Gradient uses forward differences, valid up to N-2
    EXPECT_NEAR(grad(Axis::X, 5, 5, 5), 1.0, tol);
    EXPECT_NEAR(grad(Axis::Y, 5, 5, 5), 2.0, tol);
    EXPECT_NEAR(grad(Axis::Z, 5, 5, 5), 3.0, tol);
}

TEST_F(DerivativesTest, Divergence_Composite_Test) {
    // vec = [x, 2y, 3z]
    // Divergence = 1 + 2 + 3 = 6
    // Internal points use backward difference, which is exact for linear functions.
    VectorField vec;
    vec.setup(grid,
              [](double x, double, double, double) { return x; },
              [](double, double y, double, double) { return 2.0 * y; },
              [](double, double, double z, double) { return 3.0 * z; }
    );
    vec.populate(0.0);

    Field div = createEmptyField();

    // Define BCs: Boundary values for x, 2y, 3z at their respective lower boundaries (0) are all 0.0
    Func bcFunc = [](double, double, double, double) { return 0.0; };
    double time = 0.0;

    // Updated signature: pass BC funcs for u, v, w and time
    deriv.computeDivergence(vec, div, bcFunc, bcFunc, bcFunc, time);

    // Check internal points (Backward Difference for linear is exact)
    EXPECT_NEAR(div(5, 5, 5), 6.0, 1e-8);
}

TEST_F(DerivativesTest, Laplacian_Consistency_Check) {
    // Verify: Div(Grad(f)) == Dxx(f) + Dyy(f) + Dzz(f)
    // f = x^3 + y^3 + z^3
    auto func = [](double x, double y, double z, double) {
        return x * x * x + y * y * y + z * z * z;
    };
    Field field = createField(func);

    VectorField grad = createEmptyVectorField();
    deriv.computeGradient(field, grad);

    Field divGrad = createEmptyField();

    // BCs for Divergence of Gradient:
    // This is technically evaluating the Laplacian at the boundary.
    // For consistency checking internal points, we can pass dummy BCs as long as we only check internal indices.
    Func dummyBC = [](double, double, double, double) { return 0.0; };
    double time = 0.0;

    deriv.computeDivergence(grad, divGrad, dummyBC, dummyBC, dummyBC, time);

    Field dxx = createEmptyField();
    Field dyy = createEmptyField();
    Field dzz = createEmptyField();
    deriv.computeDxx(field, dxx, dummyBC, simTime, Axis::X);
    deriv.computeDyy(field, dyy, dummyBC, simTime, Axis::Y);
    deriv.computeDzz(field, dzz, dummyBC, simTime, Axis::Z);

    for (size_t k = 1; k < Nz - 1; ++k) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            for (size_t i = 1; i < Nx - 1; ++i) {
                double sumOfSecondDerivs = dxx(i, j, k) + dyy(i, j, k) + dzz(i, j, k);
                EXPECT_NEAR(divGrad(i, j, k), sumOfSecondDerivs, 1e-10);
            }
        }
    }
}


// ============================================================================
// === EXCEPTION HANDLING
// ============================================================================

TEST_F(DerivativesTest, Throws_IfGridIsNegative) {
    EXPECT_THROW(std::make_shared<Grid>(10, 10, 10, -0.1, 1.0, 1.0), std::runtime_error);
}