#include "numerics/derivatives.hpp"

#include <stdexcept>

void Derivatives::computeGradient(const Field &field, VectorField &gradient) const {
    computeDx_fwd(field, gradient(Axis::X));
    computeDy_fwd(field, gradient(Axis::Y));
    computeDz_fwd(field, gradient(Axis::Z));
}

void Derivatives::computeDx_fwd(const Field &field, Field &dx) const {
    const auto &grid = field.getGrid();
    const double inv_dx = 1.0 / grid.dx;

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = Ny_tot * Nz_tot;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx_tot;

        // 1. Compute the derivative for each, whole x-line (last x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 0; i < Nx_tot - 1; ++i) {
            size_t idx = xLineOffset + i;
            dx[idx] = (field[idx + 1] - field[idx]) * inv_dx;
        }
    }

    // 2. Handle physical right BC: the last point of each x-line cannot compute a forward difference
    if (grid.hasMaxBoundary(Axis::X)) {
        for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
            size_t lastPhysicalXLineIdx = xLine * Nx_tot + (Nx_tot - 1 - H);
            dx[lastPhysicalXLineIdx] = 0.0;
        }
    }
}

void Derivatives::computeDy_fwd(const Field &field, Field &dy) const {
    const auto &grid = field.getGrid();
    const double inv_dy = 1.0 / grid.dy;

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // The distance in memory between (i, j, k) and (i, j+1, k)
    //  i.e., field[idx + strideY] is physically the neighbor "above" in Y w.r.t. field[idx]
    const size_t strideY = Nx_tot;
    const size_t strideZ = Ny_tot * Nx_tot;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz_tot; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (last y-elements excluded)
        for (size_t j = 0; j < Ny_tot - 1; ++j) {
            // Calculate the starting linear index of the current x-line
            const size_t xLineStart = xyPlaneOffset + (j * Nx_tot);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx_tot; ++i) {
                // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
                //  ~ "take the whole contiguous row j+1, and subtract the contiguous row j"
                size_t idx = xLineStart + i;
                dy[idx] = (field[idx + strideY] - field[idx]) * inv_dy;
            }
        }
    }

    // 2. Handle physical right BC: the last point of each y-line cannot compute a forward difference
    //  i.e., the whole last (physical) x-line of the XY-plane
    if (grid.hasMaxBoundary(Axis::Y)) {
        for (size_t k = 0; k < Nz_tot; ++k) {
            size_t xyPlaneOffset = k * strideZ;
            size_t lastPhysicalXLineStart = xyPlaneOffset + ((Ny_tot - 1 - H) * Nx_tot);
            for (size_t i = 0; i < Nx_tot; ++i) {
                dy[lastPhysicalXLineStart + i] = 0.0;
            }
        }
    }
}

void Derivatives::computeDz_fwd(const Field &field, Field &dz) const {
    const auto &grid = field.getGrid();
    const double inv_dz = 1.0 / grid.dz;

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // The distance in memory between (i, j, k) and (i, j, k+1)
    //  i.e., field[idx + strideZ] is physically the neighbor "above" in Z w.r.t. field[idx]
    const size_t strideZ = Nx_tot * Ny_tot; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, last z-elements excluded)
    for (size_t k = 0; k < Nz_tot - 1; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
            //  ~ "take the whole contiguous row z+1, and subtract the contiguous row z"
            size_t idx = xyPlaneOffset + p;
            dz[idx] = (field[idx + strideZ] - field[idx]) * inv_dz;
        }
    }

    // 2. Handle physical right BC: the last point of each z-line cannot compute a forward difference
    if (grid.hasMaxBoundary(Axis::Z)) {
        size_t lastPhysicalXyPlaneOffset = (Nz_tot - 1 - H) * strideZ;
        for (size_t p = 0; p < strideZ; ++p) {
            dz[lastPhysicalXyPlaneOffset + p] = 0.0;
        }
    }
}


void Derivatives::computeDx_bwd(const Field &field, Field &dx, Func &bcu, double time) const {
    const auto &grid = field.getGrid();
    const double inv_dx = 1.0 / grid.dx;

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Precompute constant coefficients for boundary conditions
    const double C1 = -1.0 / 3.0 * inv_dx;          // Coefficient for field(1, j, k)
    const double C2 = 3.0 * inv_dx;                 // Coefficient for field(0, j, k)
    const double C3 = -(8.0 / 3.0) * inv_dx;        // Coefficient for bcu(...)

    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = Ny_tot * Nz_tot;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx_tot;

        // 1. Compute the derivative for each, whole x-line (first x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 1; i < Nx_tot; ++i) {
            size_t idx = xLineOffset + i;
            dx[idx] = (field[idx] - field[idx - 1]) * inv_dx;
        }
    }

    // 2. Handle physical left BC: the first point of each x-line is based on the given BC-Func and time
    if (grid.hasMinBoundary(Axis::X)) {
        const size_t strideY = Nx_tot;
        const size_t strideZ = Nx_tot * Ny_tot;

        for (size_t k = 0; k < Nz_tot; ++k) {
            double physical_z = grid.to_z((long) k - (long) H, field.getOffset(), field.getOffsetAxis());
            size_t xyPlaneOffset = k * strideZ;
            for (size_t j = 0; j < Ny_tot; ++j) {
                double physical_y = grid.to_y((long) j - (long) H, field.getOffset(), field.getOffsetAxis());
                size_t firstPhysicalXLineIdx = xyPlaneOffset + j * strideY;
                double val_1 = field[firstPhysicalXLineIdx + H + 1];
                double val_0 = field[firstPhysicalXLineIdx + H];
                dx[firstPhysicalXLineIdx + H] = C1 * val_1 + C2 * val_0 + C3 * bcu(0.0, physical_y, physical_z, time);
            }
        }
    }
}

void Derivatives::computeDy_bwd(const Field &field, Field &dy, Func &bcv, double time) const {
    const auto &grid = field.getGrid();
    const double inv_dy = 1.0 / grid.dy;

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Precompute constant coefficients for boundary conditions
    const double C1 = -1.0 / 3.0 * inv_dy;          // Coefficient for field(i, 1, k)
    const double C2 = 3.0 * inv_dy;                 // Coefficient for field(i, 0, k)
    const double C3 = -(8.0 / 3.0) * inv_dy;        // Coefficient for bcv(...)

    // The distance in memory between (i, j, k) and (i, j-1, k), and (i, j, k) and (i, j, k-1)
    const size_t strideY = Nx_tot;
    const size_t strideZ = Nx_tot * Ny_tot;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz_tot; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (first y-elements excluded)
        for (size_t j = 1; j < Ny_tot; ++j) {
            const size_t xLineStart = xyPlaneOffset + (j * Nx_tot);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx_tot; ++i) {
                // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
                //  ~ "take the whole contiguous row j+1, and subtract the contiguous row j"
                size_t idx = xLineStart + i;
                dy[idx] = (field[idx] - field[idx - strideY]) * inv_dy;
            }
        }
    }

    // 2. Handle physical left BC: the first point of each y-line is based on the given BC-Func and time
    if (grid.hasMinBoundary(Axis::Y)) {
        for (size_t k = 0; k < Nz_tot; ++k) {
            double physical_z = grid.to_z((long) k - (long) H, field.getOffset(), field.getOffsetAxis());
            size_t xyPlaneOffset = k * strideZ;
            for (size_t i = 0; i < Nx_tot; ++i) {
                double physical_x = grid.to_x((long) i - (long) H, field.getOffset(), field.getOffsetAxis());
                size_t firstPhysicalXLineIdx = xyPlaneOffset + H * strideY;
                double val_1 = field[firstPhysicalXLineIdx + i + strideY];
                double val_0 = field[firstPhysicalXLineIdx + i];
                dy[firstPhysicalXLineIdx + i] = C1 * val_1 + C2 * val_0 + C3 * bcv(physical_x, 0.0, physical_z, time);
            }
        }
    }
}

void Derivatives::computeDz_bwd(const Field &field, Field &dz, Func &bcw, double time) const {
    const auto &grid = field.getGrid();
    const double inv_dz = 1.0 / grid.dz;

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Precompute constant coefficients for boundary conditions
    const double C1 = -1.0 / 3.0 * inv_dz;          // Coefficient for field(i, j, 1)
    const double C2 = 3.0 * inv_dz;                 // Coefficient for field(i, j, 0)
    const double C3 = -(8.0 / 3.0) * inv_dz;        // Coefficient for bcw(...)

    // The distance in memory between (i, j, k) and (i, j, k-1)
    const size_t strideZ = Ny_tot * Nx_tot; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, first z-elements excluded)
    for (size_t k = 1; k < Nz_tot; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 2. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
            //  ~ "take the whole contiguous row z, and subtract the contiguous row z-1"
            size_t idx = xyPlaneOffset + p;
            dz[idx] = (field[idx] - field[idx - strideZ]) * inv_dz;
        }
    }

    // 2. Handle physical left BC: the first point of each z-line (first XY-plane) is based on the given BC-Func and time
    if (grid.hasMinBoundary(Axis::Z)) {
        const size_t strideY = Nx_tot;
        const size_t firstPhysicalXyPlaneOffset = H * strideZ;
        for (size_t j = 0; j < Ny_tot; ++j) {
            const double physical_y = grid.to_y((long) j - (long) H, field.getOffset(), field.getOffsetAxis());
            const size_t xLineOffset = firstPhysicalXyPlaneOffset + (j * strideY);

            for (size_t i = 0; i < Nx_tot; ++i) {
                const double physical_x = grid.to_x((long) (i - H), field.getOffset(), field.getOffsetAxis());
                size_t idx_0 = xLineOffset + i;    // field(i, j, 0)
                size_t idx_1 = idx_0 + strideZ;    // field(i, j, 1)
                dz[idx_0] = C1 * field[idx_1] + C2 * field[idx_0] + C3 * bcw(physical_x, physical_y, 0.0, time);
            }
        }
    }
}

void Derivatives::computeDivergence(
        const VectorField &field, Field &divergence,
        Func &bcu, Func &bcv, Func &bcw, double time
) const {
    Field tmp = divergence;
    tmp.reset();
    divergence.reset();

    computeDx_bwd(field(Axis::X), tmp, bcu, time);
    divergence.add(tmp);
    computeDy_bwd(field(Axis::Y), tmp, bcv, time);
    divergence.add(tmp);
    computeDz_bwd(field(Axis::Z), tmp, bcw, time);
    divergence.add(tmp);
}


void Derivatives::computeDxx(
        const VectorField &field, VectorField &dxx,
        const Func &bcu, const Func &bcv, const Func &bcw, const double &time
) const {
    const double *fx_ptr = field(Axis::X).getData().data();
    const double *fy_ptr = field(Axis::Y).getData().data();
    const double *fz_ptr = field(Axis::Z).getData().data();
    double *dfx_dxx_ptr = dxx(Axis::X).getData().data();
    double *dfy_dxx_ptr = dxx(Axis::Y).getData().data();
    double *dfz_dxx_ptr = dxx(Axis::Z).getData().data();

    const auto &grid = field.getGrid();
    const double inv_dxx = 1.0 / (grid.dx * grid.dx);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = Ny_tot * Nz_tot;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx_tot;

        // 1. Compute the derivative for each, whole x-line (first and last x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 1; i < Nx_tot - 1; ++i) {
            size_t idx = xLineOffset + i;
            dfx_dxx_ptr[idx] = (fx_ptr[idx + 1] + fx_ptr[idx - 1] - 2.0 * fx_ptr[idx]) * inv_dxx;
            dfy_dxx_ptr[idx] = (fy_ptr[idx + 1] + fy_ptr[idx - 1] - 2.0 * fy_ptr[idx]) * inv_dxx;
            dfz_dxx_ptr[idx] = (fz_ptr[idx + 1] + fz_ptr[idx - 1] - 2.0 * fz_ptr[idx]) * inv_dxx;
        }
    }

    // 2. Handle boundary conditions
    // Precompute constant coefficients for boundary conditions
    const double C1 = 4.0 / 3.0 * inv_dxx;          // Coefficient for field(1, j, k)
    const double C2 = -4.0 * inv_dxx;               // Coefficient for field(0, j, k)
    const double C3 = (8.0 / 3.0) * inv_dxx;        // Coefficient for bc(...)

    const size_t strideY = Nx_tot;
    const size_t strideZ = Nx_tot * Ny_tot;

    // 2A. Handle physical left BC: the first point of each x-line cannot compute a centered difference
    if (grid.hasMinBoundary(Axis::X)) {
        auto &fx = field(Axis::X);
        auto fx_offset = fx.getOffset();
        auto fx_offsetAxis = fx.getOffsetAxis();

        for (size_t k = 0; k < Nz_tot; ++k) {
            double physical_Xz = grid.to_z((long) k - (long) H, fx_offset, fx_offsetAxis);
            size_t xyPlaneOffset = k * strideZ;
            for (size_t j = 0; j < Ny_tot; ++j) {
                double physical_Xy = grid.to_y((long) j - (long) H, fx_offset, fx_offsetAxis);
                size_t firstPhysicalXLineIdx = xyPlaneOffset + j * strideY;
                double val_X1 = fx_ptr[firstPhysicalXLineIdx + H + 1];
                double val_X0 = fx_ptr[firstPhysicalXLineIdx + H];
                dfx_dxx_ptr[firstPhysicalXLineIdx + H] = C1 * val_X1
                                                         + C2 * val_X0
                                                         + C3 * bcu(0.0, physical_Xy, physical_Xz, time);
                dfy_dxx_ptr[firstPhysicalXLineIdx + H] = 0;
                dfz_dxx_ptr[firstPhysicalXLineIdx + H] = 0;
            }
        }
    }
    // 2B. Handle physical right BC: the last point of each x-line cannot compute a centered difference
    if (grid.hasMaxBoundary(Axis::X)) {
        auto &fy = field(Axis::Y);
        auto fy_offset = fy.getOffset();
        auto fy_offsetAxis = fy.getOffsetAxis();

        auto &fz = field(Axis::Z);
        auto fz_offset = fz.getOffset();
        auto fz_offsetAxis = fz.getOffsetAxis();

        double physical_Xwall = grid.to_x((long)(grid.Nx - 1), GridStaggering::FACE_CENTERED, Axis::X);

        for (size_t k = 0; k < Nz_tot; ++k) {
            double physical_Yz = grid.to_z((long) k - (long) H, fy_offset, fy_offsetAxis);
            double physical_Zz = grid.to_z((long) k - (long) H, fz_offset, fz_offsetAxis);
            size_t xyPlaneOffset = k * strideZ;
            for (size_t j = 0; j < Ny_tot; ++j) {
                double physical_Yy = grid.to_y((long) j - (long) H, fy_offset, fy_offsetAxis);
                double physical_Zy = grid.to_y((long) j - (long) H, fz_offset, fz_offsetAxis);
                size_t lastPhysicalXLineIdx = (xyPlaneOffset + j * strideY) + (Nx_tot - 1 - H);
                double val_Ynm2 = fy_ptr[lastPhysicalXLineIdx - 1];
                double val_Ynm1 = fy_ptr[lastPhysicalXLineIdx];
                double val_Znm2 = fz_ptr[lastPhysicalXLineIdx - 1];
                double val_Znm1 = fz_ptr[lastPhysicalXLineIdx];
                dfx_dxx_ptr[lastPhysicalXLineIdx] = 0;
                dfy_dxx_ptr[lastPhysicalXLineIdx] = C1 * val_Ynm2
                                                    + C2 * val_Ynm1
                                                    + C3 * bcv(physical_Xwall, physical_Yy, physical_Yz, time);
                dfz_dxx_ptr[lastPhysicalXLineIdx] = C1 * val_Znm2
                                                    + C2 * val_Znm1
                                                    + C3 * bcw(physical_Xwall, physical_Zy, physical_Zz, time);
            }
        }
    }
}

void Derivatives::computeDyy(
        const VectorField &field, VectorField &dyy,
        const Func &bcu, const Func &bcv, const Func &bcw, const double &time
) const {
    const double *fx_ptr = field(Axis::X).getData().data();
    const double *fy_ptr = field(Axis::Y).getData().data();
    const double *fz_ptr = field(Axis::Z).getData().data();
    double *dfx_dyy_ptr = dyy(Axis::X).getData().data();
    double *dfy_dyy_ptr = dyy(Axis::Y).getData().data();
    double *dfz_dyy_ptr = dyy(Axis::Z).getData().data();

    const auto &grid = field.getGrid();
    const double inv_dyy = 1.0 / (grid.dy * grid.dy);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // The distance in memory between (i, j-1, k) and (i, j, k), (i, j, k) and (i, j+1, k)
    const size_t strideY = Nx_tot;
    const size_t strideZ = Ny_tot * Nx_tot;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz_tot; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (first and last y-elements excluded)
        for (size_t j = 1; j < Ny_tot - 1; ++j) {
            const size_t xLineStart = xyPlaneOffset + (j * strideY);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx_tot; ++i) {
                size_t idx = xLineStart + i;
                dfx_dyy_ptr[idx] = (fx_ptr[idx + strideY] + fx_ptr[idx - strideY] - 2.0 * fx_ptr[idx]) * inv_dyy;
                dfy_dyy_ptr[idx] = (fy_ptr[idx + strideY] + fy_ptr[idx - strideY] - 2.0 * fy_ptr[idx]) * inv_dyy;
                dfz_dyy_ptr[idx] = (fz_ptr[idx + strideY] + fz_ptr[idx - strideY] - 2.0 * fz_ptr[idx]) * inv_dyy;
            }
        }
    }

    // 2. Handle boundary conditions
    // Precompute constant coefficients for boundary conditions
    const double C1 = 4.0 / 3.0 * inv_dyy;          // Coefficient for field(i, 1, k)
    const double C2 = -4.0 * inv_dyy;               // Coefficient for field(i, 0, k)
    const double C3 = (8.0 / 3.0) * inv_dyy;        // Coefficient for bc(...)

    // 2A. Handle physical left BC: the first point of each y-line cannot compute a centered difference
    //  i.e., the whole first x-line of the XY-plane
    if (grid.hasMinBoundary(Axis::Y)) {
        auto &fy = field(Axis::Y);
        auto fy_offset = fy.getOffset();
        auto fy_offsetAxis = fy.getOffsetAxis();

        for (size_t k = 0; k < Nz_tot; ++k) {
            double physical_Yz = grid.to_z((long) k - (long) H, fy_offset, fy_offsetAxis);
            size_t xyPlaneOffset = k * strideZ;
            size_t firstPhysicalXLineIdx = xyPlaneOffset + H * strideY;
            for (size_t i = 0; i < Nx_tot; ++i) {
                double physical_Yx = grid.to_x((long) i - (long) H, fy_offset, fy_offsetAxis);
                double val_Y1 = fy_ptr[firstPhysicalXLineIdx + i + strideY];
                double val_Y0 = fy_ptr[firstPhysicalXLineIdx + i];
                dfx_dyy_ptr[firstPhysicalXLineIdx + i] = 0;
                dfy_dyy_ptr[firstPhysicalXLineIdx + i] = C1 * val_Y1
                                                         + C2 * val_Y0
                                                         + C3 * bcv(physical_Yx, 0.0, physical_Yz, time);
                dfz_dyy_ptr[firstPhysicalXLineIdx + i] = 0;
            }
        }
    }
    // 2B. Handle physical right BC: the last point of each y-line cannot compute a centered difference
    //  i.e., the whole last x-line of the XY-plane
    if (grid.hasMaxBoundary(Axis::Y)) {
        auto &fx = field(Axis::X);
        auto fx_offset = fx.getOffset();
        auto fx_offsetAxis = fx.getOffsetAxis();

        auto &fz = field(Axis::Z);
        auto fz_offset = fz.getOffset();
        auto fz_offsetAxis = fz.getOffsetAxis();

        double physical_Ywall = grid.to_y((long)(grid.Ny - 1), GridStaggering::FACE_CENTERED, Axis::Y);

        for (size_t k = 0; k < Nz_tot; ++k) {
            double physical_Xz = grid.to_z((long) k - (long) H, fx_offset, fx_offsetAxis);
            double physical_Zz = grid.to_z((long) k - (long) H, fz_offset, fz_offsetAxis);
            size_t xyPlaneOffset = k * strideZ;
            size_t lastPhysicalXLineIdx = xyPlaneOffset + (Ny_tot - 1 - H) * strideY;
            for (size_t i = 0; i < Nx_tot; ++i) {
                double physical_Xx = grid.to_x((long) i - (long) H, fx_offset, fx_offsetAxis);
                double physical_Zx = grid.to_x((long) i - (long) H, fz_offset, fz_offsetAxis);
                double val_Xnm2 = fx_ptr[lastPhysicalXLineIdx + i - strideY];
                double val_Xnm1 = fx_ptr[lastPhysicalXLineIdx + i];
                double val_Znm2 = fz_ptr[lastPhysicalXLineIdx + i - strideY];
                double val_Znm1 = fz_ptr[lastPhysicalXLineIdx + i];
                dfx_dyy_ptr[lastPhysicalXLineIdx + i] = C1 * val_Xnm2
                                                        + C2 * val_Xnm1
                                                        + C3 * bcu(physical_Xx, physical_Ywall, physical_Xz, time);
                dfy_dyy_ptr[lastPhysicalXLineIdx + i] = 0;
                dfz_dyy_ptr[lastPhysicalXLineIdx + i] = C1 * val_Znm2
                                                        + C2 * val_Znm1
                                                        + C3 * bcw(physical_Zx, physical_Ywall, physical_Zz, time);
            }
        }
    }
}

void Derivatives::computeDzz(
        const VectorField &field, VectorField &dzz,
        const Func &bcu, const Func &bcv, const Func &bcw, const double &time
) const {
    const double *fx_ptr = field(Axis::X).getData().data();
    const double *fy_ptr = field(Axis::Y).getData().data();
    const double *fz_ptr = field(Axis::Z).getData().data();
    double *dfx_dzz_ptr = dzz(Axis::X).getData().data();
    double *dfy_dzz_ptr = dzz(Axis::Y).getData().data();
    double *dfz_dzz_ptr = dzz(Axis::Z).getData().data();

    const auto &grid = field.getGrid();
    const double inv_dzz = 1.0 / (grid.dz * grid.dz);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // The distance in memory between (i, j, k-1) and (i, j, k), (i, j, k) and (i, j, k+1)
    const size_t strideZ = Nx_tot * Ny_tot; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, first and last z-elements excluded)
    for (size_t k = 1; k < Nz_tot - 1; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            size_t idx = xyPlaneOffset + p;
            dfx_dzz_ptr[idx] = (fx_ptr[idx + strideZ] + fx_ptr[idx - strideZ] - 2.0 * fx_ptr[idx]) * inv_dzz;
            dfy_dzz_ptr[idx] = (fy_ptr[idx + strideZ] + fy_ptr[idx - strideZ] - 2.0 * fy_ptr[idx]) * inv_dzz;
            dfz_dzz_ptr[idx] = (fz_ptr[idx + strideZ] + fz_ptr[idx - strideZ] - 2.0 * fz_ptr[idx]) * inv_dzz;
        }
    }

    // 2A. Handle physical left BC: the first point of each z-line cannot compute a centered difference
    if (field.getGrid().hasMinBoundary(Axis::Z)) {
        size_t firstPhysicalXyPlaneOffset = H * strideZ;
        std::fill(dfx_dzz_ptr + firstPhysicalXyPlaneOffset, dfx_dzz_ptr + firstPhysicalXyPlaneOffset + strideZ, 0.0);
        std::fill(dfy_dzz_ptr + firstPhysicalXyPlaneOffset, dfy_dzz_ptr + firstPhysicalXyPlaneOffset + strideZ, 0.0);
        std::fill(dfz_dzz_ptr + firstPhysicalXyPlaneOffset, dfz_dzz_ptr + firstPhysicalXyPlaneOffset + strideZ, 0.0);
    }
    // 2B. Handle physical right BC: the last point of each z-line cannot compute a centered difference
    if (field.getGrid().hasMaxBoundary(Axis::Z)) {
        size_t lastPhysicalXyPlaneOffset = (Nz_tot - 1 - H) * strideZ;
        std::fill(dfx_dzz_ptr + lastPhysicalXyPlaneOffset, dfx_dzz_ptr + lastPhysicalXyPlaneOffset + strideZ, 0.0);
        std::fill(dfy_dzz_ptr + lastPhysicalXyPlaneOffset, dfy_dzz_ptr + lastPhysicalXyPlaneOffset + strideZ, 0.0);
        std::fill(dfz_dzz_ptr + lastPhysicalXyPlaneOffset, dfz_dzz_ptr + lastPhysicalXyPlaneOffset + strideZ, 0.0);
    }
}

void Derivatives::computeDxx(
        const Field &field, Field &dxx,
        const Func &bc, const double &time, const Axis axis
) const {
    const auto &grid = field.getGrid();
    const double inv_dxx = 1.0 / (grid.dx * grid.dx);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = Ny_tot * Nz_tot;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx_tot;

        // 1. Compute the derivative for each, whole x-line (first and last x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 1; i < Nx_tot - 1; ++i) {
            size_t idx = xLineOffset + i;
            dxx[idx] = (field[idx + 1] + field[idx - 1] - 2.0 * field[idx]) * inv_dxx;
        }
    }

    // 2. Handle boundary conditions
    // Precompute constant coefficients for boundary conditions
    const double C1 = 4.0 / 3.0 * inv_dxx;          // Coefficient for field(1, j, k)
    const double C2 = -4.0 * inv_dxx;               // Coefficient for field(0, j, k)
    const double C3 = (8.0 / 3.0) * inv_dxx;        // Coefficient for bc(...)

    const size_t strideY = Nx_tot;
    const size_t strideZ = Nx_tot * Ny_tot;

    auto offset = field.getOffset();
    auto offsetAxis = field.getOffsetAxis();

    if (axis == Axis::X) {
        // 2A. Handle physical left BC: the first point of each x-line cannot compute a centered difference
        if (grid.hasMinBoundary(Axis::X)) {
            for (size_t k = 0; k < Nz_tot; ++k) {
                double physical_z = grid.to_z((long) k - (long) H, offset, offsetAxis);
                size_t xyPlaneOffset = k * strideZ;
                for (size_t j = 0; j < Ny_tot; ++j) {
                    double physical_y = grid.to_y((long) j - (long) H, offset, offsetAxis);
                    size_t firstPhysicalXLineIdx = xyPlaneOffset + j * strideY;
                    double val_1 = field[firstPhysicalXLineIdx + H + 1];
                    double val_0 = field[firstPhysicalXLineIdx + H];
                    dxx[firstPhysicalXLineIdx + H] = C1 * val_1
                                                     + C2 * val_0
                                                     + C3 * bc(0.0, physical_y, physical_z, time);
                }
            }
        }
        // 2B. Handle physical right BC: the last point of each x-line cannot compute a centered difference
        if (grid.hasMaxBoundary(Axis::X)) {
            for (size_t k = 0; k < Nz_tot; ++k) {
                size_t xyPlaneOffset = k * strideZ;
                for (size_t j = 0; j < Ny_tot; ++j) {
                    size_t lastPhysicalXLineIdx = (xyPlaneOffset + j * strideY) + (Nx_tot - 1 - H);
                    dxx[lastPhysicalXLineIdx] = 0.0;
                }
            }
        }
    } else { // not an X-component
        // 2A. Handle physical left BC: the first point of each x-line cannot compute a centered difference
        if (grid.hasMinBoundary(Axis::X)) {
            for (size_t k = 0; k < Nz_tot; ++k) {
                size_t xyPlaneOffset = k * strideZ;
                for (size_t j = 0; j < Ny_tot; ++j) {
                    size_t firstPhysicalXLineIdx = xyPlaneOffset + j * strideY;
                    dxx[firstPhysicalXLineIdx + H] = 0;
                }
            }
        }
        // 2B. Handle physical right BC: the last point of each x-line cannot compute a centered difference
        if (grid.hasMaxBoundary(Axis::X)) {
            double physical_Xwall = grid.to_x((long)(grid.Nx - 1), GridStaggering::FACE_CENTERED, Axis::X);
            for (size_t k = 0; k < Nz_tot; ++k) {
                double physical_z = grid.to_z((long) k - (long) H, offset, offsetAxis);
                size_t xyPlaneOffset = k * strideZ;
                for (size_t j = 0; j < Ny_tot; ++j) {
                    double physical_y = grid.to_y((long) j - (long) H, offset, offsetAxis);
                    size_t lastPhysicalXLineIdx = (xyPlaneOffset + j * strideY) + (Nx_tot - 1 - H);
                    double val_nm2 = field[lastPhysicalXLineIdx - 1];
                    double val_nm1 = field[lastPhysicalXLineIdx];
                    dxx[lastPhysicalXLineIdx] = C1 * val_nm2
                                                + C2 * val_nm1
                                                + C3 * bc(physical_Xwall, physical_y, physical_z, time);
                }
            }
        }
    }
}

void Derivatives::computeDyy(
        const Field &field, Field &dyy,
        const Func &bc, const double &time, const Axis axis
) const {
    const auto &grid = field.getGrid();
    const double inv_dyy = 1.0 / (grid.dy * grid.dy);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // The distance in memory between (i, j-1, k) and (i, j, k), (i, j, k) and (i, j+1, k)
    const size_t strideY = Nx_tot;
    const size_t strideZ = Ny_tot * Nx_tot;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz_tot; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (first and last y-elements excluded)
        for (size_t j = 1; j < Ny_tot - 1; ++j) {
            const size_t xLineStart = xyPlaneOffset + (j * strideY);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx_tot; ++i) {
                size_t idx = xLineStart + i;
                dyy[idx] = (field[idx + strideY] + field[idx - strideY] - 2.0 * field[idx]) * inv_dyy;
            }
        }
    }

    // 2. Handle boundary conditions
    // Precompute constant coefficients for boundary conditions
    const double C1 = 4.0 / 3.0 * inv_dyy;          // Coefficient for field(i, 1, k)
    const double C2 = -4.0 * inv_dyy;               // Coefficient for field(i, 0, k)
    const double C3 = (8.0 / 3.0) * inv_dyy;        // Coefficient for bc(...)

    auto offset = field.getOffset();
    auto offsetAxis = field.getOffsetAxis();

    if (axis == Axis::Y) {
        // 2A. Handle physical left BC: the first point of each y-line cannot compute a centered difference
        //  i.e., the whole first x-line of the XY-plane
        if (grid.hasMinBoundary(Axis::Y)) {
            for (size_t k = 0; k < Nz_tot; ++k) {
                double physical_z = grid.to_z((long) k - (long) H, offset, offsetAxis);
                size_t xyPlaneOffset = k * strideZ;
                size_t firstPhysicalXLineIdx = xyPlaneOffset + H * strideY;
                for (size_t i = 0; i < Nx_tot; ++i) {
                    double physical_x = grid.to_x((long) i - (long) H, offset, offsetAxis);
                    double val_1 = field[firstPhysicalXLineIdx + i + strideY];
                    double val_0 = field[firstPhysicalXLineIdx + i];
                    dyy[firstPhysicalXLineIdx + i] = C1 * val_1
                                                     + C2 * val_0
                                                     + C3 * bc(physical_x, 0.0, physical_z, time);
                }
            }
        }
        // 2B. Handle physical right BC: the last point of each y-line cannot compute a centered difference
        //  i.e., the whole last x-line of the XY-plane
        if (grid.hasMaxBoundary(Axis::Y)) {
            for (size_t k = 0; k < Nz_tot; ++k) {
                size_t xyPlaneOffset = k * strideZ;
                size_t lastPhysicalXLineIdx = xyPlaneOffset + (Ny_tot - 1 - H) * strideY;
                for (size_t i = 0; i < Nx_tot; ++i) {
                    dyy[lastPhysicalXLineIdx + i] = 0;
                }
            }
        }
    } else {
        // 2A. Handle physical left BC: the first point of each y-line cannot compute a centered difference
        //  i.e., the whole first x-line of the XY-plane
        if (grid.hasMinBoundary(Axis::Y)) {
            for (size_t k = 0; k < Nz_tot; ++k) {
                size_t xyPlaneOffset = k * strideZ;
                size_t firstPhysicalXLineIdx = xyPlaneOffset + H * strideY;
                for (size_t i = 0; i < Nx_tot; ++i) {
                    dyy[firstPhysicalXLineIdx + i] = 0;
                }
            }
        }
        // 2B. Handle physical right BC: the last point of each y-line cannot compute a centered difference
        //  i.e., the whole last x-line of the XY-plane
        if (grid.hasMaxBoundary(Axis::Y)) {
            double physical_Ywall = grid.to_y((long)(grid.Ny - 1), GridStaggering::FACE_CENTERED, Axis::Y);
            for (size_t k = 0; k < Nz_tot; ++k) {
                double physical_z = grid.to_z((long) k - (long) H, offset, offsetAxis);
                size_t xyPlaneOffset = k * strideZ;
                size_t lastPhysicalXLineIdx = xyPlaneOffset + (Ny_tot - 1 - H) * strideY;
                for (size_t i = 0; i < Nx_tot; ++i) {
                    double physical_x = grid.to_x((long) i - (long) H, offset, offsetAxis);
                    double val_nm2 = field[lastPhysicalXLineIdx + i - strideY];
                    double val_nm1 = field[lastPhysicalXLineIdx + i];
                    dyy[lastPhysicalXLineIdx + i] = C1 * val_nm2
                                                    + C2 * val_nm1
                                                    + C3 * bc(physical_x, physical_Ywall, physical_z, time);
                }
            }
        }
    }
}

void Derivatives::computeDzz(
        const Field &field, Field &dzz,
        const Func &bc, const double &time, const Axis axis
) const {
    const auto &grid = field.getGrid();
    const double inv_dzz = 1.0 / (grid.dz * grid.dz);

    const size_t H = grid.n_halo;
    const size_t Nx_tot = grid.Nx + 2 * H;
    const size_t Ny_tot = grid.Ny + 2 * H;
    const size_t Nz_tot = grid.Nz + 2 * H;

    // The distance in memory between (i, j, k-1) and (i, j, k), (i, j, k) and (i, j, k+1)
    const size_t strideZ = Nx_tot * Ny_tot; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, first and last z-elements excluded)
    for (size_t k = 1; k < Nz_tot - 1; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            size_t idx = xyPlaneOffset + p;
            dzz[idx] = (field[idx + strideZ] + field[idx - strideZ] - 2.0 * field[idx]) * inv_dzz;
        }
    }

    // 2A. Handle physical left BC: the first point of each z-line cannot compute a centered difference
    if (grid.hasMinBoundary(Axis::Z)) {
        for (size_t p = 0; p < strideZ; ++p) {
            size_t firstPhysicalXyPlaneIdx = H * strideZ + p;
            dzz[firstPhysicalXyPlaneIdx] = 0.0;
        }
    }
    // 2B. Handle physical right BC: the last point of each z-line cannot compute a centered difference
    if (grid.hasMaxBoundary(Axis::Z)) {
        for (size_t p = 0; p < strideZ; ++p) {
            size_t lastPhysicalXyPlaneIdx = (Nz_tot - 1 - H) * strideZ + p;
            dzz[lastPhysicalXyPlaneIdx] = 0.0;
        }
    }

}
