#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>

#include "core/Fields.hpp"
#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"

TEST(DerivativesTest, ComputeDx_QuadraticField_CorrectDerivatives) {
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
