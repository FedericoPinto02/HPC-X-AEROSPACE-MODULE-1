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
    const size_t lastPlaneOffset = (Nz - 1) * strideZ;
    for (size_t p = 0; p < strideZ; ++p) {
        dz[lastPlaneOffset + p] = 0.0;
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

void Derivatives::computeHessianDiag(const Field &field, VectorField &hessianDiag) const {
    computeDxx(field, hessianDiag(Axis::X));
    computeDyy(field, hessianDiag(Axis::Y));
    computeDzz(field, hessianDiag(Axis::Z));
}

void Derivatives::computeDxx(const Field &field, Field &dxx) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / (grid.dx * grid.dx);

    if (grid.Nx < 3) {
        throw std::runtime_error("grid size Nx must be at least 3 to compute second derivative.");
    }

    for (size_t k = 0; k < grid.Nz; k++) {
        for (size_t j = 0; j < grid.Ny; j++) {
            for (size_t i = 1; i < (grid.Nx - 1); i++) {
                dxx(i, j, k) = (field(i + 1, j, k) + field(i - 1, j, k) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }
}

void Derivatives::computeDyy(const Field &field, Field &dyy) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / (grid.dy * grid.dy);

    if (grid.Ny < 3) {
        throw std::runtime_error("grid size Ny must be at least 3 to compute second derivative.");
    }

    for (size_t k = 0; k < grid.Nz; k++) {
        for (size_t j = 1; j < (grid.Ny - 1); j++) {
            for (size_t i = 0; i < grid.Nx; i++) {
                dyy(i, j, k) = (field(i, j + 1, k) + field(i, j - 1, k) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }
}

void Derivatives::computeDzz(const Field &field, Field &dzz) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / (grid.dz * grid.dz);

    if (grid.Nz < 3) {
        throw std::runtime_error("grid size Nz must be at least 3 to compute second derivative.");
    }

    for (size_t k = 1; k < (grid.Nz - 1); k++) {
        for (size_t j = 0; j < grid.Ny; j++) {
            for (size_t i = 0; i < grid.Nx; i++) {
                dzz(i, j, k) = (field(i, j, k + 1) + field(i, j, k - 1) - 2.0 * field(i, j, k)) * mul;
            }
        }
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