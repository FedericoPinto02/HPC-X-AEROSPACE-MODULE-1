#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>

#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"

TEST(ComputeDx, ComputeDx_QuadraticField_CorrectDerivatives) {
    // === Grid setup ===
    std::shared_ptr<const Grid> grid = std::make_shared<Grid>(10, 1, 1, 0.1, 1.0, 1.0);

    // === Field: f(x) = x^2 ===
    std::vector<Field::Scalar> values(grid->size());
    size_t idx = 0;
    for (size_t k = 0; k < grid->Nz; ++k)
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t i = 0; i < grid->Nx; ++i, ++idx)
                values[idx] = (i * grid->dx) * (i * grid->dx); // x^2

    // === Setup Field ===
    Field field;
    field.setup(grid, values);

    // === Setup Field for derivative ===
    Field dfdx;
    dfdx.setup(grid, std::vector<Field::Scalar>(grid->size(), 0.0));

    // === Compute derivative ===
    Derivatives deriv;
    deriv.computeDx(field, dfdx);

    // === Check derivatives in interior ===
    const double tol = 1e-6;
    for (size_t i = 1; i < grid->Nx - 1; ++i) {
        double x = i * grid->dx;
        double expected = 2.0 * x; // analytic derivative
        EXPECT_NEAR(dfdx(i, 0, 0), expected, tol) << "Error at i=" << i;
    }

    // === Borders should stay at 0 ===
    EXPECT_DOUBLE_EQ(dfdx(0, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(dfdx(grid->Nx - 1, 0, 0), 0.0);
}

TEST(ComputeDx, Throws_NullGrid) {
    Field field;
    Field dx;

    // field and dx not set up → grid = nullptr
    Derivatives deriv;

    EXPECT_THROW(deriv.computeDx(field, dx), std::runtime_error);
}

TEST(ComputeDx, Throws_NxDim){
    std::shared_ptr<const Grid> grid = std::make_shared<Grid>(2,1,1,0.1,1,1);
    Field field, dx;
    field.setup(grid, std::vector<double>(grid->size(), 1.0));
    dx.setup(grid, std::vector<double>(grid->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDx(field, dx), std::runtime_error);
}

TEST(ComputeDx, Throws_FieldGridDiff) {
    std::shared_ptr<const Grid> grid1 = std::make_shared<Grid>(10,1,1,0.1,1,1);
    std::shared_ptr<const Grid> grid2 = std::make_shared<Grid>(10,1,1,0.1,1,1);

    Field field;
    field.setup(grid1, std::vector<double>(grid1->size(), 1.0));

    Field dx;
    dx.setup(grid2, std::vector<double>(grid2->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDx(field, dx), std::runtime_error);
}

TEST(ComputeDx, Throws_NegativeDx) {
    std::shared_ptr<const Grid> grid = std::make_shared<Grid>(10,1,1,-0.1,1,1); // dx < 0 
    Field field, dx;
    field.setup(grid, std::vector<double>(grid->size(), 1.0));
    dx.setup(grid, std::vector<double>(grid->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDx(field, dx), std::runtime_error);
}

TEST(ComputeDy, ComputeDy_MixedFunctionXY2_CorrectDerivatives) {
    // === Grid setup ===
    std::shared_ptr<const Grid> grid = std::make_shared<Grid>(5, 5, 1, 0.5, 0.5, 1.0);

    // === Field: f(x,y) = x * y^2 ===
    std::vector<Field::Scalar> values(grid->size());
    size_t idx = 0;
    for (size_t k = 0; k < grid->Nz; ++k)
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t i = 0; i < grid->Nx; ++i, ++idx) {
                double x = i * grid->dx;
                double y = j * grid->dy;
                values[idx] = x * y * y;
            }

    // === Setup Field and derivative ===
    Field field;
    field.setup(grid, values);

    Field dfdy;
    dfdy.setup(grid, std::vector<Field::Scalar>(grid->size(), 0.0));

    // === Compute derivative ===
    Derivatives deriv;
    deriv.computeDy(field, dfdy);

    // === Check derivatives in interior ===
    const double tol = 1e-6;
    for (size_t j = 1; j < grid->Ny - 1; ++j)
        for (size_t i = 0; i < grid->Nx; ++i) {
            double x = i * grid->dx;
            double y = j * grid->dy;
            double expected = 2.0 * x * y;
            EXPECT_NEAR(dfdy(i, j, 0), expected, tol)
                << "Error at (i,j)=(" << i << "," << j << ")";
        }
    
    // === Borders should stay at 0 ===
    EXPECT_DOUBLE_EQ(dfdy(0, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(dfdy(0, grid->Ny - 1, 0), 0.0);
}

TEST(ComputeDy, ComputeDy_Throws_NullGrid) {
    Field field;
    Field dy;

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDy(field, dy), std::runtime_error);
}

TEST(ComputeDy, ComputeDy_Throws_FieldGridDiff) {
    auto grid1 = std::make_shared<Grid>(1,10,1,1.0,0.1,1.0);
    auto grid2 = std::make_shared<Grid>(1,10,1,1.0,0.1,1.0);

    Field field;
    field.setup(grid1, std::vector<double>(grid1->size(), 1.0));

    Field dy;
    dy.setup(grid2, std::vector<double>(grid2->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDy(field, dy), std::runtime_error);
}

TEST(ComputeDy, ComputeDy_Throws_NyTooSmall) {
    auto grid = std::make_shared<Grid>(1,2,1,1.0,0.1,1.0); // Ny=2
    Field field, dy;
    field.setup(grid, std::vector<double>(grid->size(), 1.0));
    dy.setup(grid, std::vector<double>(grid->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDy(field, dy), std::runtime_error);
}

TEST(ComputeDy, ComputeDy_Throws_NegativeDy) {
    auto grid = std::make_shared<Grid>(1,10,1,1.0,-0.1,1.0); // dy negativo
    Field field, dy;
    field.setup(grid, std::vector<double>(grid->size(), 1.0));
    dy.setup(grid, std::vector<double>(grid->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDy(field, dy), std::runtime_error);
}

