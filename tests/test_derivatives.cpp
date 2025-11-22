#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <functional>

#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"

class DerivativesTest : public ::testing::Test {
protected:
    // Constants
    const size_t Nx = 10, Ny = 10, Nz = 10;
    const double dx = 0.5, dy = 0.5, dz = 0.5;

    // Stack-allocated Grid
    const Grid grid;
    Derivatives deriv;

    DerivativesTest() : grid(Nx, Ny, Nz, dx, dy, dz) {}

    // Helper to quickly setup and populate a field with a lambda
    Field createField(Functions::Func func) {
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
// === CONSISTENCY TEST (Integration)
// ============================================================================

TEST_F(DerivativesTest, Laplacian_Consistency) {
    // Ideally: Div(Grad(f)) == Laplacian(f)
    // In your implementation:
    // Grad is Forward Difference. Div is Backward Difference.
    // Forward + Backward == Central Difference, which is exactly your Hessian implementation.

    // f(x,y,z) = x^3 + y^3 + z^3
    auto func = [](double x, double y, double z, double) {
        return x*x*x + y*y*y + z*z*z;
    };
    Field field = createField(func);

    // 1. Compute Div(Grad(f))
    VectorField grad = createEmptyVectorField();
    deriv.computeGradient(field, grad);

    Field divGrad = createEmptyField();
    deriv.computeDivergence(grad, divGrad);

    // 2. Compute Hessian Diagonals directly (Central Diff)
    VectorField hessian = createEmptyVectorField();
    deriv.computeHessianDiag(field, hessian);

    // 3. Compare Sum of Hessian Diags vs Div(Grad)
    // Valid only in the interior [1, N-2] because boundaries handle 0s differently
    const double tol = 1e-10;
    for (size_t k = 1; k < Nz - 1; ++k) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            for (size_t i = 1; i < Nx - 1; ++i) {
                double laplacianVal = hessian(Axis::X, i, j, k) +
                                      hessian(Axis::Y, i, j, k) +
                                      hessian(Axis::Z, i, j, k);

                EXPECT_NEAR(divGrad(i, j, k), laplacianVal, tol)
                                    << "Mismatch at " << i << "," << j << "," << k;
            }
        }
    }
}

// ============================================================================
// === FIRST DERIVATIVE TESTS (Forward Difference)
// ============================================================================

TEST_F(DerivativesTest, Gradient_ForwardDifference_Accuracy) {
    // f(x,y,z) = x^2 + y^2 + z^2
    Field field = createField([](double x, double y, double z, double) {
        return x*x + y*y + z*z;
    });

    VectorField grad = createEmptyVectorField();
    deriv.computeGradient(field, grad);

    // Forward Difference of x^2 is ( (x+h)^2 - x^2 ) / h = 2x + h
    const double tol = 1e-6;
    for (size_t k = 0; k < Nz - 1; ++k) {
        for (size_t j = 0; j < Ny - 1; ++j) {
            for (size_t i = 0; i < Nx - 1; ++i) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;

                EXPECT_NEAR(grad(Axis::X, i, j, k), 2.0 * x + dx, tol);
                EXPECT_NEAR(grad(Axis::Y, i, j, k), 2.0 * y + dy, tol);
                EXPECT_NEAR(grad(Axis::Z, i, j, k), 2.0 * z + dz, tol);
            }
        }
    }
}

TEST_F(DerivativesTest, Gradient_Boundaries_AreZero) {
    Field field = createField(Functions::ZERO); // Content doesn't matter for boundary check
    VectorField grad = createEmptyVectorField();
    deriv.computeGradient(field, grad);

    // Verify last index is explicitly set to 0.0 as per implementation
    for (size_t k = 0; k < Nz; ++k)
        for (size_t j = 0; j < Ny; ++j)
            EXPECT_DOUBLE_EQ(grad(Axis::X, Nx - 1, j, k), 0.0);
}

// ============================================================================
// === DIVERGENCE TESTS (Backward Difference)
// ============================================================================

TEST_F(DerivativesTest, Divergence_BackwardDifference_Accuracy) {
    // u = (x, 2y, 3z)
    VectorField vec = createEmptyVectorField();
    vec.component(Axis::X).setup(grid, [](double x, double, double, double){ return x; });
    vec.component(Axis::Y).setup(grid, [](double, double y, double, double){ return 2.0*y; });
    vec.component(Axis::Z).setup(grid, [](double, double, double z, double){ return 3.0*z; });
    vec.populate(0.0);

    Field div = createEmptyField();
    deriv.computeDivergence(vec, div);

    // Backward diff of linear function (ax) is exactly 'a'.
    // Div = 1 + 2 + 3 = 6.
    // Valid for indices [1, N-1]
    const double tol = 1e-6;
    for (size_t k = 1; k < Nz; ++k) {
        for (size_t j = 1; j < Ny; ++j) {
            for (size_t i = 1; i < Nx; ++i) {
                EXPECT_NEAR(div(i, j, k), 6.0, tol);
            }
        }
    }
}

TEST_F(DerivativesTest, Divergence_Boundaries_AreZero) {
    VectorField vec = createEmptyVectorField();
    Field div = createEmptyField();
    deriv.computeDivergence(vec, div);

    // Verify first index is explicitly set to 0.0 as per implementation
    for (size_t k = 0; k < Nz; ++k)
        for (size_t j = 0; j < Ny; ++j)
            EXPECT_DOUBLE_EQ(div(0, j, k), 0.0);
}

// ============================================================================
// === SECOND DERIVATIVE TESTS (Central Difference)
// ============================================================================

TEST_F(DerivativesTest, Hessian_CentralDifference_ExactForCubic) {
    // f(x,y,z) = x^3 + y^3 + z^3
    Field field = createField([](double x, double y, double z, double) {
        return x*x*x + y*y*y + z*z*z;
    });

    VectorField hessian = createEmptyVectorField();
    deriv.computeHessianDiag(field, hessian);

    // Central Difference is exact for cubic polynomials.
    // d2/dx2 (x^3) = 6x
    const double tol = 1e-8;
    for (size_t k = 1; k < Nz - 1; ++k) {
        for (size_t j = 1; j < Ny - 1; ++j) {
            for (size_t i = 1; i < Nx - 1; ++i) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;

                EXPECT_NEAR(hessian(Axis::X, i, j, k), 6.0 * x, tol);
                EXPECT_NEAR(hessian(Axis::Y, i, j, k), 6.0 * y, tol);
                EXPECT_NEAR(hessian(Axis::Z, i, j, k), 6.0 * z, tol);
            }
        }
    }
}

// ============================================================================
// === EXCEPTION HANDLING
// ============================================================================

TEST_F(DerivativesTest, Throws_IfGridIsTooSmall) {
    Grid tinyGrid(2, 2, 2, 1.0, 1.0, 1.0); // Too small for 2nd derivative (needs 3 points)
    Field field;
    field.setup(tinyGrid);
    Field result;
    result.setup(tinyGrid);

    EXPECT_THROW(deriv.computeDxx(field, result), std::runtime_error);
    EXPECT_THROW(deriv.computeDyy(field, result), std::runtime_error);
    EXPECT_THROW(deriv.computeDzz(field, result), std::runtime_error);
}

TEST_F(DerivativesTest, Throws_IfGridIsNegative) {
    EXPECT_THROW(Grid(10, 10, 10, -0.1, 1.0, 1.0), std::runtime_error);
}