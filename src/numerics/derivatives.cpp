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
    const double mul = 1.0 / grid.dx;

    for (size_t k = 0; k < grid.Nz; k++) {
        for (size_t j = 0; j < grid.Ny; j++) {
            for (size_t i = 0; i < (grid.Nx - 1); i++) {
                dx(i, j, k) = (field(i + 1, j, k) - field(i, j, k)) * mul;
            }
            dx(grid.Nx - 1, j, k) = 0.0;
        }
    }
}

void Derivatives::computeDy_fwd(const Field &field, Field &dy) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dy;

    for (size_t k = 0; k < grid.Nz; k++) {
        for (size_t i = 0; i < grid.Nx; i++) {
            for (size_t j = 0; j < (grid.Ny - 1); j++) {
                dy(i, j, k) = (field(i, j + 1, k) - field(i, j, k)) * mul;
            }
            dy(i, grid.Ny - 1, k) = 0.0;
        }
    }
}

void Derivatives::computeDz_fwd(const Field &field, Field &dz) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dz;

    for (size_t i = 0; i < grid.Nx; i++) {
        for (size_t j = 0; j < grid.Ny; j++) {
            for (size_t k = 0; k < (grid.Nz - 1); k++) {
                dz(i, j, k) = (field(i, j, k + 1) - field(i, j, k)) * mul;
            }
            dz(i, j, grid.Nz - 1) = 0.0;
        }
    }
}


void Derivatives::computeDx_bwd(const Field &field, Field &dx) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dx;

    for (size_t k = 0; k < grid.Nz; k++) {
        for (size_t j = 0; j < grid.Ny; j++) {
            for (size_t i = 1; i < grid.Nx; i++) {
                dx(i, j, k) = (field(i, j, k) - field(i - 1, j, k)) * mul;
            }
            dx(0, j, k) = 0.0;
        }
    }
}

void Derivatives::computeDy_bwd(const Field &field, Field &dy) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dy;

    for (size_t k = 0; k < grid.Nz; k++) {
        for (size_t i = 0; i < grid.Nx; i++) {
            for (size_t j = 1; j < grid.Ny; j++) {
                dy(i, j, k) = (field(i, j, k) - field(i, j - 1, k)) * mul;
            }
            dy(i, 0, k) = 0.0;
        }
    }
}

void Derivatives::computeDz_bwd(const Field &field, Field &dz) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dz;

    for (size_t i = 0; i < grid.Nx; i++) {
        for (size_t j = 0; j < grid.Ny; j++) {
            for (size_t k = 1; k < grid.Nz; k++) {
                dz(i, j, k) = (field(i, j, k) - field(i, j, k - 1)) * mul;
            }
            dz(i, j, 0) = 0.0;
        }
    }
}


void Derivatives::computeDivergence(const VectorField &field, Field &divergence) const {
    Field tmp = divergence;
    tmp.reset();
    divergence.reset();

    computeDx_bwd(field(Axis::X), tmp);
    divergence.add(tmp);
    computeDy_bwd(field(Axis::Y), tmp);
    divergence.add(tmp);
    computeDz_bwd(field(Axis::Z), tmp);
    divergence.add(tmp);
}

void Derivatives::computeHessianDiag(const Field &field, VectorField &hessianDiag) const {
    computeDxx(field, hessianDiag(Axis::X));
    computeDyy(field, hessianDiag(Axis::Y));
    computeDzz(field, hessianDiag(Axis::Z));
}

void Derivatives::computeDxx(const VectorField &field, VectorField &dxx) const {
    computeDxx(field(Axis::X), dxx(Axis::X));
    computeDxx(field(Axis::Y), dxx(Axis::Y));
    computeDxx(field(Axis::Z), dxx(Axis::Z));
}

void Derivatives::computeDyy(const VectorField &field, VectorField &dyy) const {
    computeDyy(field(Axis::X), dyy(Axis::X));
    computeDyy(field(Axis::Y), dyy(Axis::Y));
    computeDyy(field(Axis::Z), dyy(Axis::Z));
}

void Derivatives::computeDzz(const VectorField &field, VectorField &dzz) const {
    computeDzz(field(Axis::X), dzz(Axis::X));
    computeDzz(field(Axis::Y), dzz(Axis::Y));
    computeDzz(field(Axis::Z), dzz(Axis::Z));
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