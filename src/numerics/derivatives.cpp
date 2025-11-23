#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"
#include <stdexcept>

void Derivatives::computeGradient(const Field &field, VectorField &gradient) const {
    computeDx(field, gradient(Axis::X));
    computeDy(field, gradient(Axis::Y));
    computeDz(field, gradient(Axis::Z));
}

void Derivatives::computeDx(const Field &field, Field &dx) const {
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

void Derivatives::computeDy(const Field &field, Field &dy) const {
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

void Derivatives::computeDz(const Field &field, Field &dz) const {
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


void Derivatives::computeDxDiv(const Field &field, Field &dx) const {
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

void Derivatives::computeDyDiv(const Field &field, Field &dy) const {
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

void Derivatives::computeDzDiv(const Field &field, Field &dz) const {
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

    computeDxDiv(field(Axis::X), tmp);
    divergence.add(tmp);
    computeDyDiv(field(Axis::Y), tmp);
    divergence.add(tmp);
    computeDzDiv(field(Axis::Z), tmp);
    divergence.add(tmp);
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