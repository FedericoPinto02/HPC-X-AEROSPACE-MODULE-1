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

TEST(ComputeDz, ComputeDz_MixedFunctionXYZ_CorrectDerivatives) {
    // === Grid setup ===
    std::shared_ptr<const Grid> grid = std::make_shared<Grid>(3, 3, 5, 1.0, 1.0, 0.5);

    // === Field: f(x,y,z) = x * y * z^2 ===
    std::vector<Field::Scalar> values(grid->size());
    size_t idx = 0;
    for (size_t k = 0; k < grid->Nz; ++k)
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t i = 0; i < grid->Nx; ++i, ++idx) {
                double x = i * grid->dx;
                double y = j * grid->dy;
                double z = k * grid->dz;
                values[idx] = x * y * z * z;
            }

    // === Setup Field and derivative ===
    Field field;
    field.setup(grid, values);

    Field dfdz;
    dfdz.setup(grid, std::vector<Field::Scalar>(grid->size(), 0.0));

    // === Compute derivative ===
    Derivatives deriv;
    deriv.computeDz(field, dfdz);

    // === Check derivatives in interior ===
    const double tol = 1e-6;
    for (size_t k = 1; k < grid->Nz - 1; ++k)
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t i = 0; i < grid->Nx; ++i) {
                double x = i * grid->dx;
                double y = j * grid->dy;
                double z = k * grid->dz;
                double expected = 2.0 * x * y * z;
                EXPECT_NEAR(dfdz(i, j, k), expected, tol)
                    << "Error at (i,j,k)=(" << i << "," << j << "," << k << ")";
            }

    // === Borders should stay at 0 ===
    for (size_t j = 0; j < grid->Ny; ++j)
        for (size_t i = 0; i < grid->Nx; ++i) {
            EXPECT_DOUBLE_EQ(dfdz(i, j, 0), 0.0);
            EXPECT_DOUBLE_EQ(dfdz(i, j, grid->Nz - 1), 0.0);
        }
}

TEST(ComputeDz, ComputeDz_Throws_NullGrid) {
    Field field;
    Field dz;

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDz(field, dz), std::runtime_error);
}

TEST(ComputeDz, ComputeDz_Throws_FieldGridDiff) {
    auto grid1 = std::make_shared<Grid>(3,3,5,1.0,1.0,0.5);
    auto grid2 = std::make_shared<Grid>(3,3,5,1.0,1.0,0.5);

    Field field;
    field.setup(grid1, std::vector<double>(grid1->size(), 1.0));

    Field dz;
    dz.setup(grid2, std::vector<double>(grid2->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDz(field, dz), std::runtime_error);
}

TEST(ComputeDz, ComputeDz_Throws_NzTooSmall) {
    auto grid = std::make_shared<Grid>(3,3,2,1.0,1.0,0.5); // Nz=2
    Field field, dz;
    field.setup(grid, std::vector<double>(grid->size(), 1.0));
    dz.setup(grid, std::vector<double>(grid->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDz(field, dz), std::runtime_error);
}

TEST(ComputeDz, ComputeDz_Throws_NegativeDz) {
    auto grid = std::make_shared<Grid>(3,3,5,1.0,1.0,-0.5); // dz negativo
    Field field, dz;
    field.setup(grid, std::vector<double>(grid->size(), 1.0));
    dz.setup(grid, std::vector<double>(grid->size(), 0.0));

    Derivatives deriv;
    EXPECT_THROW(deriv.computeDz(field, dz), std::runtime_error);
}

TEST(GeneralDerivativeTest, ComputeAllDerivatives) {
    // === Grid setup ===
    std::shared_ptr<const Grid> grid = std::make_shared<Grid>(5, 5, 5, 0.5, 0.5, 0.5);

    // === Field: f(x,y,z) = x^2 + y^2 + z^2 ===
    std::vector<Field::Scalar> values(grid->size());
    size_t idx = 0;
    for (size_t k = 0; k < grid->Nz; ++k)
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t i = 0; i < grid->Nx; ++i, ++idx) {
                double x = i * grid->dx;
                double y = j * grid->dy;
                double z = k * grid->dz;
                values[idx] = x*x + y*y + z*z;
            }

    // === Setup Field and derivatives ===
    Field field;
    field.setup(grid, values);

    Field dfdx, dfdy, dfdz;
    dfdx.setup(grid, std::vector<double>(grid->size(), 0.0));
    dfdy.setup(grid, std::vector<double>(grid->size(), 0.0));
    dfdz.setup(grid, std::vector<double>(grid->size(), 0.0));

    // === Compute derivatives ===
    Derivatives deriv;
    deriv.computeDx(field, dfdx);
    deriv.computeDy(field, dfdy);
    deriv.computeDz(field, dfdz);

    // === Check derivatives in interior ===
    const double tol = 1e-6;
    for (size_t k = 1; k < grid->Nz - 1; ++k)
        for (size_t j = 1; j < grid->Ny - 1; ++j)
            for (size_t i = 1; i < grid->Nx - 1; ++i) {
                double x = i * grid->dx;
                double y = j * grid->dy;
                double z = k * grid->dz;

                EXPECT_NEAR(dfdx(i,j,k), 2.0*x, tol)
                    << "Error in dfdx at (i,j,k)=(" << i << "," << j << "," << k << ")";
                EXPECT_NEAR(dfdy(i,j,k), 2.0*y, tol)
                    << "Error in dfdy at (i,j,k)=(" << i << "," << j << "," << k << ")";
                EXPECT_NEAR(dfdz(i,j,k), 2.0*z, tol)
                    << "Error in dfdz at (i,j,k)=(" << i << "," << j << "," << k << ")";
            }

    // === Borders should stay at 0 ===
    for (size_t k : {size_t(0), grid->Nz-1})
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t i = 0; i < grid->Nx; ++i) {
                EXPECT_DOUBLE_EQ(dfdz(i,j,k), 0.0);
            }

    for (size_t j : {size_t(0), grid->Ny-1})
        for (size_t k = 0; k < grid->Nz; ++k)
            for (size_t i = 0; i < grid->Nx; ++i) {
                EXPECT_DOUBLE_EQ(dfdy(i,j,k), 0.0);
            }

    for (size_t i : {size_t(0), grid->Nx-1})
        for (size_t j = 0; j < grid->Ny; ++j)
            for (size_t k = 0; k < grid->Nz; ++k) {
                EXPECT_DOUBLE_EQ(dfdx(i,j,k), 0.0);
            }
}
