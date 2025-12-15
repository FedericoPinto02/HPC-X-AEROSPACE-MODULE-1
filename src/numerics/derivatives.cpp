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

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                double val_i = field(i, j, k);
                double val_ip1 = field.valueWithOffset(i, j, k, Axis::X, +1);
                dx(i, j, k) = (val_ip1 - val_i) * mul;
            }
        }
    }

    // Handle physical right boundary: zero
    if (grid.hasMaxBoundary(Axis::X)) {
        long i = (long) grid.Nx - 1;
        for (long k = 0; k < grid.Nz; k++) {
            for (long j = 0; j < grid.Ny; j++) {
                dx(i, j, k) = 0.0;
            }
        }
    }
}

void Derivatives::computeDy_fwd(const Field &field, Field &dy) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dy;

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                double val_j = field(i, j, k);
                double val_jp1 = field.valueWithOffset(i, j, k, Axis::Y, +1);
                dy(i, j, k) = (val_jp1 - val_j) * mul;
            }
        }
    }

    // Handle physical right boundary: zero
    if (grid.hasMaxBoundary(Axis::Y)) {
        long j = grid.Ny - 1;
        for (long k = 0; k < grid.Nz; k++) {
            for (long i = 0; i < grid.Nx; i++) {
                dy(i, j, k) = 0.0;
            }
        }
    }
}

void Derivatives::computeDz_fwd(const Field &field, Field &dz) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dz;

    for (long i = 0; i < grid.Nx; i++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long k = 0; k < grid.Nz; k++) {
                double val_k = field(i, j, k);
                double val_kp1 = field.valueWithOffset(i, j, k, Axis::Z, +1);
                dz(i, j, k) = (val_kp1 - val_k) * mul;
            }
        }
    }

    // Handle physical right boundary: zero
    if (grid.hasMaxBoundary(Axis::Z)) {
        long k = (long) grid.Nz - 1;
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dz(i, j, k) = 0.0;
            }
        }
    }
}


void Derivatives::computeDx_bwd(const Field &field, Field &dx, Func &bcu, double time) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dx;

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dx(i, j, k) = (field(i, j, k) - field(i - 1, j, k)) * mul;
            }
        }
    }

    // Handle physical left boundary: Taylor expansion with ghost point
    if (grid.hasMinBoundary(Axis::X)) {
        for (long k = 0; k < grid.Nz; k++) {
            double physical_z = grid.to_z(k, field.getOffset(), field.getOffsetAxis());
            for (long j = 0; j < grid.Ny; j++) {
                double physical_y = grid.to_y(j, field.getOffset(), field.getOffsetAxis());
                double val_1 = field(1, j, k);
                double val_0 = field(0, j, k);
                dx(0, j, k) = -1.0 / 3.0 * mul * val_1 + 3.0 * mul * val_0
                              - (8.0 / 3.0) * mul * bcu(0.0, physical_y, physical_z, time);
            }
        }
    }
}

void Derivatives::computeDy_bwd(const Field &field, Field &dy, Func &bcv, double time) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dy;

    for (long k = 0; k < grid.Nz; k++) {
        for (long i = 0; i < grid.Nx; i++) {
            for (long j = 0; j < grid.Ny; j++) {
                dy(i, j, k) = (field(i, j, k) - field(i, j - 1, k)) * mul;
            }
        }
    }

    // Handle physical left boundary: Taylor expansion with ghost point
    if (grid.hasMinBoundary(Axis::Y)) {
        for (long k = 0; k < grid.Nz; k++) {
            double physical_z = grid.to_z(k, field.getOffset(), field.getOffsetAxis());
            for (long i = 0; i < grid.Nx; i++) {
                double physical_x = grid.to_x(i, field.getOffset(), field.getOffsetAxis());
                double val_1 = field(i, 1, k);
                double val_0 = field(i, 0, k);
                dy(i, 0, k) = -1.0 / 3.0 * mul * val_1 + 3.0 * mul * val_0
                              - (8.0 / 3.0) * mul * bcv(physical_x, 0.0, physical_z, time);
            }
        }
    }
}

void Derivatives::computeDz_bwd(const Field &field, Field &dz, Func &bcw, double time) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / grid.dz;

    for (long i = 0; i < grid.Nx; i++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long k = 0; k < grid.Nz; k++) {
                dz(i, j, k) = (field(i, j, k) - field(i, j, k - 1)) * mul;
            }
        }
    }

    // Handle physical left boundary: Taylor expansion with ghost point
    if (grid.hasMinBoundary(Axis::Z)) {
        for (long j = 0; j < grid.Ny; j++) {
            double physical_y = grid.to_y(j, field.getOffset(), field.getOffsetAxis());
            for (long i = 0; i < grid.Nx; i++) {
                double physical_x = grid.to_x(i, field.getOffset(), field.getOffsetAxis());
                double val_1 = field(i, j, 1);
                double val_0 = field(i, j, 0);
                dz(i, j, 0) = -1.0 / 3.0 * mul * val_1 + 3.0 * mul * val_0
                              - (8.0 / 3.0) * mul * bcw(physical_x, physical_y, 0.0, time);
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

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dxx(i, j, k) = (field(i + 1, j, k) + field(i - 1, j, k) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }

    if (grid.hasMinBoundary(Axis::X)) {
        for (long k = 0; k < grid.Nz; k++) {
            for (long j = 0; j < grid.Ny; j++) {
                dxx(0, j, k) = 0.0;
            }
        }
    }
    if (grid.hasMaxBoundary(Axis::X)) {
        long i = (long) grid.Nx - 1;
        for (long k = 0; k < grid.Nz; k++) {
            for (long j = 0; j < grid.Ny; j++) {
                dxx(i, j, k) = 0.0;
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

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dyy(i, j, k) = (field(i, j + 1, k) + field(i, j - 1, k) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }

    if (grid.hasMinBoundary(Axis::Y)) {
        for (long k = 0; k < grid.Nz; k++) {
            for (long i = 0; i < grid.Nx; i++) {
                dyy(i, 0, k) = 0.0;
            }
        }
    }
    if (grid.hasMaxBoundary(Axis::Y)) {
        long j = (long) grid.Ny - 1;
        for (long k = 0; k < grid.Nz; k++) {
            for (long i = 0; i < grid.Nx; i++) {
                dyy(i, j, k) = 0.0;
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

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dzz(i, j, k) = (field(i, j, k + 1) + field(i, j, k - 1) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }

    if (grid.hasMinBoundary(Axis::Z)) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dzz(i, j, 0) = 0.0;
            }
        }
    }
    if (grid.hasMaxBoundary(Axis::Z)) {
        long k = (long) grid.Nz - 1;
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dzz(i, j, k) = 0.0;
            }
        }
    }
}

double Derivatives::Dxx_local(const Field &f, size_t i, size_t j, size_t k) const {
    const auto &grid = f.getGrid();
    const double mul = 1.0 / (grid.dx * grid.dx);

    if (i == 0 || i == grid.Nx - 1) return 0.0; // o BC

    return (f((long) i + 1, (long) j, (long) k)
            + f((long) i - 1, (long) j, (long) k)
            - 2.0 * f((long) i, (long) j, (long) k)) * mul;
}

double Derivatives::Dyy_local(const Field &f, size_t i, size_t j, size_t k) const {
    const auto &grid = f.getGrid();
    const double mul = 1.0 / (grid.dy * grid.dy);

    if (j == 0 || j == grid.Ny - 1) return 0.0; // o BC

    return (f((long) i, (long) j + 1, (long) k)
            + f((long) i, (long) j - 1, (long) k)
            - 2.0 * f((long) i, (long) j, (long) k)) * mul;
}

double Derivatives::Dzz_local(const Field &f, size_t i, size_t j, size_t k) const {
    const auto &grid = f.getGrid();
    const double mul = 1.0 / (grid.dz * grid.dz);

    if (k == 0 || k == grid.Nz - 1) return 0.0; // o BC
    return (f((long) i, (long) j, (long) k + 1)
            + f((long) i, (long) j, (long) k - 1)
            - 2.0 * f((long) i, (long) j, (long) k)) * mul;
}