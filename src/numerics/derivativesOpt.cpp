#include "numerics/derivatives.hpp"
#include "core/Grid.hpp"
#include <stdexcept>

void Derivatives::computeGradient(const Field &field, VectorField &gradient) const {
    computeDx_fwd(field, gradient(Axis::X));
    computeDy_fwd(field, gradient(Axis::Y));
    computeDz_fwd(field, gradient(Axis::Z));
}

void Derivatives::computeDx_fwd(const Field &field, Field &dx) const {
    const auto &grid = field.getGrid();
    const double inv_dx = 1.0 / grid.dx;

    const size_t Nx = grid.Nx;
    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = grid.Ny * grid.Nz;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx;

        // 1. Compute the derivative for each, whole x-line (last x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 0; i < Nx - 1; ++i) {
            size_t idx = xLineOffset + i;
            dx[idx] = (field[idx + 1] - field[idx]) * inv_dx;
        }

        // 2. Handle Boundary Condition: the last point of the x-line cannot compute a forward difference
        dx[xLineOffset + Nx - 1] = 0.0;
    }
}

void Derivatives::computeDy_fwd(const Field &field, Field &dy) const {
    const auto &grid = field.getGrid();
    const double inv_dy = 1.0 / grid.dy;

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;

    // The distance in memory between (i, j, k) and (i, j+1, k)
    //  i.e., field[idx + strideY] is physically the neighbor "above" in Y w.r.t. field[idx]
    const size_t strideY = Nx;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz; ++k) {
        const size_t xyPlaneOffset = k * Ny * Nx; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (last y-elements excluded)
        for (size_t j = 0; j < Ny - 1; ++j) {
            // Calculate the starting linear index of the current x-line
            const size_t xLineStart = xyPlaneOffset + (j * Nx);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx; ++i) {
                // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
                //  ~ "take the whole contiguous row j+1, and subtract the contiguous row j"
                size_t idx = xLineStart + i;
                dy[idx] = (field[idx + strideY] - field[idx]) * inv_dy;
            }
        }

        // 2. Handle Boundary Condition: the last point of each y-line cannot compute a forward difference
        //  i.e., the whole last x-line of the XY-plane
        const size_t lastXLineStart = xyPlaneOffset + ((Ny - 1) * Nx);
        for (size_t i = 0; i < Nx; ++i) {
            dy[lastXLineStart + i] = 0.0;
        }
    }
}

void Derivatives::computeDz_fwd(const Field &field, Field &dz) const {
    const auto &grid = field.getGrid();
    const double inv_dz = 1.0 / grid.dz;

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;

    // The distance in memory between (i, j, k) and (i, j, k+1)
    //  i.e., field[idx + strideZ] is physically the neighbor "above" in Z w.r.t. field[idx]
    const size_t strideZ = Nx * Ny; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, last z-elements excluded)
    for (size_t k = 0; k < Nz - 1; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
            //  ~ "take the whole contiguous row z+1, and subtract the contiguous row z"
            size_t idx = xyPlaneOffset + p;
            dz[idx] = (field[idx + strideZ] - field[idx]) * inv_dz;
        }
    }

    // 2. Handle Boundary Condition: the last point of each z-line cannot compute a forward difference
    const size_t lastXyPlaneOffset = (Nz - 1) * strideZ;
    for (size_t p = 0; p < strideZ; ++p) {
        dz[lastXyPlaneOffset + p] = 0.0;
    }
}


void Derivatives::computeDx_bwd(const Field &field, Field &dx) const {
    const auto &grid = field.getGrid();
    const double inv_dx = 1.0 / grid.dx;

    const size_t Nx = grid.Nx;
    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = grid.Ny * grid.Nz;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx;

        // 1. Compute the derivative for each, whole x-line (first x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 1; i < Nx; ++i) {
            size_t idx = xLineOffset + i;
            dx[idx] = (field[idx] - field[idx - 1]) * inv_dx;
        }

        // 2. Handle Boundary Condition: the first point of the x-line cannot compute a backward difference
        dx[xLineOffset] = 0.0;
    }
}

void Derivatives::computeDy_bwd(const Field &field, Field &dy) const {
    const auto &grid = field.getGrid();
    const double inv_dy = 1.0 / grid.dy;

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;

    // The distance in memory between (i, j, k) and (i, j-1, k)
    //  i.e., field[idx - strideY] is physically the neighbor "below" in Y w.r.t. field[idx]
    const size_t strideY = Nx;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz; ++k) {
        const size_t xyPlaneOffset = k * Ny * Nx; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (first y-elements excluded)
        for (size_t j = 1; j < Ny; ++j) {
            // Calculate the starting linear index of the current x-line
            const size_t xLineStart = xyPlaneOffset + (j * Nx);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx; ++i) {
                // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
                //  ~ "take the whole contiguous row j, and subtract the contiguous row j-1"
                size_t idx = xLineStart + i;
                dy[idx] = (field[idx] - field[idx - strideY]) * inv_dy;
            }
        }

        // 2. Handle Boundary Condition: the first point of each y-line cannot compute a backward difference
        //  i.e., the whole first x-line of the XY-plane
        for (size_t i = 0; i < Nx; ++i) {
            dy[xyPlaneOffset + i] = 0.0;
        }
    }
}

void Derivatives::computeDz_bwd(const Field &field, Field &dz) const {
    const auto &grid = field.getGrid();
    const double inv_dz = 1.0 / grid.dz;

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;

    // The distance in memory between (i, j, k) and (i, j, k-1)
    //  i.e., field[idx - strideZ] is physically the neighbor "below" in Z w.r.t. field[idx]
    const size_t strideZ = Nx * Ny; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, first z-elements excluded)
    for (size_t k = 1; k < Nz; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            // Striding ensures we are reading/writing contiguous memory blocks in an efficient way
            //  ~ "take the whole contiguous row z, and subtract the contiguous row z-1"
            size_t idx = xyPlaneOffset + p;
            dz[idx] = (field[idx] - field[idx - strideZ]) * inv_dz;
        }
    }

    // 2. Handle Boundary Condition: the first point of each z-line cannot compute a backward difference
    for (size_t p = 0; p < strideZ; ++p) {
        dz[p] = 0.0;
    }
}


void Derivatives::computeDivergence(const VectorField &field, Field &divergence) const {
    Field tmp = divergence; // to clone shared-pointer grid
    divergence.reset();

    computeDx_bwd(field(Axis::X), tmp);
    divergence.add(tmp);
    computeDy_bwd(field(Axis::Y), tmp);
    divergence.add(tmp);
    computeDz_bwd(field(Axis::Z), tmp);
    divergence.add(tmp);

    // todo - could be further optimized to avoid time-interleaved addition, reset and temporary field allocation
}


void Derivatives::computeDxx(const VectorField &field, VectorField &dxx) const {
    const double *fx_ptr = field(Axis::X).getData().data();
    const double *fy_ptr = field(Axis::Y).getData().data();
    const double *fz_ptr = field(Axis::Z).getData().data();
    double *dfx_dxx_ptr = dxx(Axis::X).getData().data();
    double *dfy_dxx_ptr = dxx(Axis::Y).getData().data();
    double *dfz_dxx_ptr = dxx(Axis::Z).getData().data();

    const auto &grid = field.getGrid();
    const double inv_dxx = 1.0 / (grid.dx * grid.dx);

    const size_t Nx = grid.Nx;
    if (grid.Ny < 3) {
        throw std::runtime_error("Grid size Ny must be at least 3 to compute second derivative.");
    }

    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = grid.Ny * grid.Nz;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx;

        // 1. Compute the derivative for each, whole x-line (first and last x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 1; i < Nx - 1; ++i) {
            size_t idx = xLineOffset + i;
            dfx_dxx_ptr[idx] = (fx_ptr[idx + 1] + fx_ptr[idx - 1] - 2.0 * fx_ptr[idx]) * inv_dxx;
            dfy_dxx_ptr[idx] = (fy_ptr[idx + 1] + fy_ptr[idx - 1] - 2.0 * fy_ptr[idx]) * inv_dxx;
            dfz_dxx_ptr[idx] = (fz_ptr[idx + 1] + fz_ptr[idx - 1] - 2.0 * fz_ptr[idx]) * inv_dxx;
        }

        // 2. Handle Boundary Condition: the first and last points of the x-line cannot compute a centered difference
        dfx_dxx_ptr[xLineOffset] = 0.0;
        dfx_dxx_ptr[xLineOffset + Nx - 1] = 0.0;
        dfy_dxx_ptr[xLineOffset] = 0.0;
        dfy_dxx_ptr[xLineOffset + Nx - 1] = 0.0;
        dfz_dxx_ptr[xLineOffset] = 0.0;
        dfz_dxx_ptr[xLineOffset + Nx - 1] = 0.0;
    }
}

void Derivatives::computeDyy(const VectorField &field, VectorField &dyy) const {
    const double *fx_ptr = field(Axis::X).getData().data();
    const double *fy_ptr = field(Axis::Y).getData().data();
    const double *fz_ptr = field(Axis::Z).getData().data();
    double *dfx_dyy_ptr = dyy(Axis::X).getData().data();
    double *dfy_dyy_ptr = dyy(Axis::Y).getData().data();
    double *dfz_dyy_ptr = dyy(Axis::Z).getData().data();

    const auto &grid = field.getGrid();
    const double inv_dyy = 1.0 / (grid.dy * grid.dy);

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;
    if (grid.Ny < 3) {
        throw std::runtime_error("Grid size Nz must be at least 3 to compute second derivative.");
    }

    // The distance in memory between (i, j-1, k) and (i, j, k), (i, j, k) and (i, j+1, k)
    const size_t strideY = Nx;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz; ++k) {
        const size_t xyPlaneOffset = k * Ny * Nx; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (first and last y-elements excluded)
        for (size_t j = 1; j < Ny - 1; ++j) {
            const size_t xLineStart = xyPlaneOffset + (j * Nx);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx; ++i) {
                size_t idx = xLineStart + i;
                dfx_dyy_ptr[idx] = (fx_ptr[idx + strideY] + fx_ptr[idx - strideY] - 2.0 * fx_ptr[idx]) * inv_dyy;
                dfy_dyy_ptr[idx] = (fy_ptr[idx + strideY] + fy_ptr[idx - strideY] - 2.0 * fy_ptr[idx]) * inv_dyy;
                dfz_dyy_ptr[idx] = (fz_ptr[idx + strideY] + fz_ptr[idx - strideY] - 2.0 * fz_ptr[idx]) * inv_dyy;
            }
        }

        // 2. Handle Boundary Condition: the first and last points of each y-line cannot compute a centered difference
        //  i.e., the whole first and last x-lines of the XY-plane
        const size_t lastXLineStart = xyPlaneOffset + ((Ny - 1) * Nx);
        std::fill(dfx_dyy_ptr + xyPlaneOffset, dfx_dyy_ptr + xyPlaneOffset + Nx, 0.0);      // first x-line
        std::fill(dfx_dyy_ptr + lastXLineStart, dfx_dyy_ptr + lastXLineStart + Nx, 0.0);    // last x-line
        std::fill(dfy_dyy_ptr + xyPlaneOffset, dfy_dyy_ptr + xyPlaneOffset + Nx, 0.0);      // first x-line
        std::fill(dfy_dyy_ptr + lastXLineStart, dfy_dyy_ptr + lastXLineStart + Nx, 0.0);    // last x-line
        std::fill(dfz_dyy_ptr + xyPlaneOffset, dfz_dyy_ptr + xyPlaneOffset + Nx, 0.0);      // first x-line
        std::fill(dfz_dyy_ptr + lastXLineStart, dfz_dyy_ptr + lastXLineStart + Nx, 0.0);    // last x-line
    }
}

void Derivatives::computeDzz(const VectorField &field, VectorField &dzz) const {
    const double *fx_ptr = field(Axis::X).getData().data();
    const double *fy_ptr = field(Axis::Y).getData().data();
    const double *fz_ptr = field(Axis::Z).getData().data();
    double *dfx_dzz_ptr = dzz(Axis::X).getData().data();
    double *dfy_dzz_ptr = dzz(Axis::Y).getData().data();
    double *dfz_dzz_ptr = dzz(Axis::Z).getData().data();

    const auto &grid = field.getGrid();
    const double inv_dzz = 1.0 / (grid.dz * grid.dz);

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;
    if (grid.Nz < 3) {
        throw std::runtime_error("Grid size Nz must be at least 3 to compute second derivative.");
    }

    // The distance in memory between (i, j, k-1) and (i, j, k), (i, j, k) and (i, j, k+1)
    const size_t strideZ = Nx * Ny; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, first and last z-elements excluded)
    for (size_t k = 1; k < Nz - 1; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            size_t idx = xyPlaneOffset + p;
            dfx_dzz_ptr[idx] = (fx_ptr[idx + strideZ] + fx_ptr[idx - strideZ] - 2.0 * fx_ptr[idx]) * inv_dzz;
            dfy_dzz_ptr[idx] = (fy_ptr[idx + strideZ] + fy_ptr[idx - strideZ] - 2.0 * fy_ptr[idx]) * inv_dzz;
            dfz_dzz_ptr[idx] = (fz_ptr[idx + strideZ] + fz_ptr[idx - strideZ] - 2.0 * fz_ptr[idx]) * inv_dzz;
        }
    }

    // 2. Handle Boundary Condition: the first and last points of each z-line cannot compute a centered difference
    const size_t lastXyPlaneOffset = (Nz - 1) * strideZ;
    std::fill(dfx_dzz_ptr, dfx_dzz_ptr + strideZ, 0.0);                                         // first XY-plane
    std::fill(dfx_dzz_ptr + lastXyPlaneOffset, dfx_dzz_ptr + lastXyPlaneOffset + strideZ, 0.0); // last XY-plane
    std::fill(dfy_dzz_ptr, dfy_dzz_ptr + strideZ, 0.0);                                         // first XY-plane
    std::fill(dfy_dzz_ptr + lastXyPlaneOffset, dfy_dzz_ptr + lastXyPlaneOffset + strideZ, 0.0); // last XY-plane
    std::fill(dfz_dzz_ptr, dfz_dzz_ptr + strideZ, 0.0);                                         // first XY-plane
    std::fill(dfz_dzz_ptr + lastXyPlaneOffset, dfz_dzz_ptr + lastXyPlaneOffset + strideZ, 0.0); // last XY-plane
}

void Derivatives::computeDxx(const Field &field, Field &dxx) const {
    const auto &grid = field.getGrid();
    const double inv_dxx = 1.0 / (grid.dx * grid.dx);

    const size_t Nx = grid.Nx;
    if (grid.Ny < 3) {
        throw std::runtime_error("Grid size Ny must be at least 3 to compute second derivative.");
    }

    // Flatten Y and Z dimensions into a single 'xLine' count
    const size_t totalXLines = grid.Ny * grid.Nz;

    // Iterate over every x-line in the volume
    for (size_t xLine = 0; xLine < totalXLines; ++xLine) {
        const size_t xLineOffset = xLine * Nx;

        // 1. Compute the derivative for each, whole x-line (first and last x-elements excluded)
        // Linear index iteration: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
        for (size_t i = 1; i < Nx - 1; ++i) {
            size_t idx = xLineOffset + i;
            dxx[idx] = (field[idx + 1] + field[idx - 1] - 2.0 * field[idx]) * inv_dxx;
        }

        // 2. Handle Boundary Condition: the first and last points of the x-line cannot compute a centered difference
        dxx[xLineOffset] = 0.0;
        dxx[xLineOffset + Nx - 1] = 0.0;
    }
}

void Derivatives::computeDyy(const Field &field, Field &dyy) const {
    const auto &grid = field.getGrid();
    const double inv_dyy = 1.0 / (grid.dy * grid.dy);

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;
    if (grid.Ny < 3) {
        throw std::runtime_error("Grid size Nz must be at least 3 to compute second derivative.");
    }

    // The distance in memory between (i, j-1, k) and (i, j, k), (i, j, k) and (i, j+1, k)
    const size_t strideY = Nx;

    // Iterate over all XY-planes (at each vertical level k)
    for (size_t k = 0; k < Nz; ++k) {
        const size_t xyPlaneOffset = k * Ny * Nx; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (first and last y-elements excluded)
        for (size_t j = 1; j < Ny - 1; ++j) {
            const size_t xLineStart = xyPlaneOffset + (j * Nx);

            // Still iterate over each x-line, with linear indexes: maximum cache coherency, allowing compiler auto-vectorization (SIMD)
            for (size_t i = 0; i < Nx; ++i) {
                size_t idx = xLineStart + i;
                dyy[idx] = (field[idx + strideY] + field[idx - strideY] - 2.0 * field[idx]) * inv_dyy;
            }
        }

        // 2. Handle Boundary Condition: the first and last points of each y-line cannot compute a centered difference
        //  i.e., the whole first and last x-lines of the XY-plane
        const size_t lastXLineStart = xyPlaneOffset + ((Ny - 1) * Nx);
        for (size_t i = 0; i < Nx; ++i) {
            dyy[xyPlaneOffset + i] = 0.0;    // first x-line
            dyy[lastXLineStart + i] = 0.0;   // last x-line
        }
    }
}

void Derivatives::computeDzz(const Field &field, Field &dzz) const {
    const auto &grid = field.getGrid();
    const double inv_dzz = 1.0 / (grid.dz * grid.dz);

    const size_t Nx = grid.Nx;
    const size_t Ny = grid.Ny;
    const size_t Nz = grid.Nz;
    if (grid.Nz < 3) {
        throw std::runtime_error("Grid size Nz must be at least 3 to compute second derivative.");
    }

    // The distance in memory between (i, j, k-1) and (i, j, k), (i, j, k) and (i, j, k+1)
    const size_t strideZ = Nx * Ny; // stride is one full XY plane

    // Iterate over all XY-planes (at each vertical level k, first and last z-elements excluded)
    for (size_t k = 1; k < Nz - 1; ++k) {
        const size_t xyPlaneOffset = k * strideZ; // skip an entire XY plane for each k-level

        // 1. Compute the derivative for each, whole XY-plane (flattened: each point, one by one)
        for (size_t p = 0; p < strideZ; ++p) {
            size_t idx = xyPlaneOffset + p;
            dzz[idx] = (field[idx + strideZ] + field[idx - strideZ] - 2.0 * field[idx]) * inv_dzz;
        }
    }

    // 2. Handle Boundary Condition: the first and last points of each z-line cannot compute a centered difference
    const size_t lastXyPlaneOffset = (Nz - 1) * strideZ;
    for (size_t p = 0; p < strideZ; ++p) {
        dzz[p] = 0.0;                     // first XY-plane
        dzz[lastXyPlaneOffset + p] = 0.0; // last XY-plane
    }
}


double Derivatives::Dxx_local(const Field &f, size_t i, size_t j, size_t k) const {
    const auto &grid = f.getGrid();
    const double mul = 1.0 / (grid.dx * grid.dx);

    if (i == 0 || i == grid.Nx - 1) return 0.0; // o BC

    return (f(i + 1, j, k) + f(i - 1, j, k) - 2.0 * f(i, j, k)) * mul;
}

double Derivatives::Dyy_local(const Field &f, size_t i, size_t j, size_t k) const {
    const auto &grid = f.getGrid();
    const double mul = 1.0 / (grid.dy * grid.dy);

    if (j == 0 || j == grid.Ny - 1) return 0.0; // o BC

    return (f(i, j + 1, k) + f(i, j - 1, k) - 2.0 * f(i, j, k)) * mul;
}

double Derivatives::Dzz_local(const Field &f, size_t i, size_t j, size_t k) const {
    const auto &grid = f.getGrid();
    const double mul = 1.0 / (grid.dz * grid.dz);

    if (k == 0 || k == grid.Nz - 1) return 0.0; // o BC
    return (f(i, j, k + 1) + f(i, j, k - 1) - 2.0 * f(i, j, k)) * mul;
}