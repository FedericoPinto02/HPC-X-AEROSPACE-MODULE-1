#include "numerics/derivatives.hpp"

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


void Derivatives::computeDxx(const VectorField &field, VectorField &dxx, const Func &bcu, const Func &bcv, const Func &bcw, const double &time) const {
    computeDxx(field(Axis::X), dxx(Axis::X), bcu, time, Axis::X);
    computeDxx(field(Axis::Y), dxx(Axis::Y), bcv, time, Axis::Y);
    computeDxx(field(Axis::Z), dxx(Axis::Z), bcw, time, Axis::Z);
}

void Derivatives::computeDyy(const VectorField &field, VectorField &dyy, const Func &bcu, const Func &bcv, const Func &bcw, const double &time) const {
    computeDyy(field(Axis::X), dyy(Axis::X), bcu, time, Axis::X);
    computeDyy(field(Axis::Y), dyy(Axis::Y), bcv, time, Axis::Y);
    computeDyy(field(Axis::Z), dyy(Axis::Z), bcw, time, Axis::Z);
}

void Derivatives::computeDzz(const VectorField &field, VectorField &dzz, const Func &bcu, const Func &bcv, const Func &bcw, const double &time) const {
    computeDzz(field(Axis::X), dzz(Axis::X), bcu, time, Axis::X);
    computeDzz(field(Axis::Y), dzz(Axis::Y), bcv, time, Axis::Y);
    computeDzz(field(Axis::Z), dzz(Axis::Z), bcw, time, Axis::Z);
}

void Derivatives::computeDxx(const Field &field, Field &dxx, const Func &bc, const double &time, const Axis axis) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / (grid.dx * grid.dx);


    for (long k = 0; k < grid.Nz; k++)
    {
        for (long j = 0; j < grid.Ny; j++)
        {
            for (long i = 0; i < grid.Nx; i++)
            {
                dxx(i, j, k) = (field(i + 1, j, k) + field(i - 1, j, k) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }
    
    if (axis == Axis::X) {
        if (grid.hasMinBoundary(Axis::X)) {
            for (long k = 0; k < grid.Nz; k++) {
                double physical_z = grid.to_z(k, field.getOffset(), field.getOffsetAxis());
                for (long j = 0; j < grid.Ny; j++) {
                    double physical_y = grid.to_y(j, field.getOffset(), field.getOffsetAxis());
                    dxx(0, j, k) =( 4.0 / 3.0 *  field(1, j, k) - 4.0 * field(0, j, k) 
                                 + (8.0 / 3.0) *  bc(0.0, physical_y, physical_z, time)) * mul;
                }
            }
        }

        if (grid.hasMaxBoundary(Axis::X)) {
            for (long k = 0; k < grid.Nz; k++) {
                for (long j = 0; j < grid.Ny; j++) {
                    dxx(grid.Nx - 1, j, k) = 0.0;
                }
            }
        }
    } else {
        if (grid.hasMinBoundary(Axis::X)) {
            for (long k = 0; k < grid.Nz; k++) {
                for (long j = 0; j < grid.Ny; j++) {
                    dxx(0, j, k) = 0.0;
                }
            }
        }

        if (grid.hasMaxBoundary(Axis::X)) {
            for (long k = 0; k < grid.Nz; k++) {
                double physical_z = grid.to_z(k, field.getOffset(), field.getOffsetAxis());
                for (long j = 0; j < grid.Ny; j++) {
                    double physical_y = grid.to_y(j, field.getOffset(), field.getOffsetAxis());
                    dxx(grid.Nx - 1, j, k) = (4.0 / 3.0 * field(grid.Nx - 2, j, k) 
                                           - 4.0 *  field(grid.Nx - 1, j, k) 
                                           + (8.0 / 3.0) *  bc((grid.Nx - 0.5) * grid.dx, physical_y, physical_z, time)) * mul;
                }
            }
        }
    }
}

void Derivatives::computeDyy(const Field &field, Field &dyy, const Func &bc, const double &time, const Axis axis) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / (grid.dy * grid.dy);

    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dyy(i, j, k) = (field(i, j + 1, k) + field(i, j - 1, k) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }

    if (axis == Axis::Y) {
        if (grid.hasMinBoundary(Axis::Y)) {
            for (long k = 0; k < grid.Nz; k++) {
                double physical_z = grid.to_z(k, field.getOffset(), field.getOffsetAxis());
                for (long i = 0; i < grid.Nx; i++) {
                    double physical_x = grid.to_x(i, field.getOffset(), field.getOffsetAxis());
                    dyy(i, 0, k) = (4.0 / 3.0 *  field(i, 1, k) - 4.0 *  field(i, 0, k) 
                                 + (8.0 / 3.0) *  bc(physical_x, 0.0, physical_z, time)) * mul;
                }
            }
        }

        if (grid.hasMaxBoundary(Axis::Y)) {
            for (long k = 0; k < grid.Nz; k++) {
                for (long i = 0; i < grid.Nx; i++) {
                    dyy(i, grid.Ny - 1, k) = 0.0;
                }
            }
        }
    } else {
        if (grid.hasMinBoundary(Axis::Y)) {
            for (long k = 0; k < grid.Nz; k++) {
                for (long i = 0; i < grid.Nx; i++) {
                    dyy(i, 0, k) = 0.0;
                }
            }
        }

        if (grid.hasMaxBoundary(Axis::Y)) {
            for (long k = 0; k < grid.Nz; k++) {
                double physical_z = grid.to_z(k, field.getOffset(), field.getOffsetAxis());
                for (long i = 0; i < grid.Nx; i++) {
                    double physical_x = grid.to_x(i, field.getOffset(), field.getOffsetAxis());
                    dyy(i, grid.Ny - 1, k) = (4.0 / 3.0 *  field(i, grid.Ny - 2, k) 
                                           - 4.0 *  field(i, grid.Ny - 1, k) 
                                           + (8.0 / 3.0) *  bc(physical_x, (grid.Ny - 0.5) * grid.dy, physical_z, time)) * mul;
                }
            }
        }
    }
}

void Derivatives::computeDzz(const Field &field, Field &dzz, const Func &bc, const double &time, const Axis axis) const {
    const auto &grid = field.getGrid();
    const double mul = 1.0 / (grid.dz * grid.dz);


    for (long k = 0; k < grid.Nz; k++) {
        for (long j = 0; j < grid.Ny; j++) {
            for (long i = 0; i < grid.Nx; i++) {
                dzz(i, j, k) = (field(i, j, k + 1) + field(i, j, k - 1) - 2.0 * field(i, j, k)) * mul;
            }
        }
    }

    if (axis == Axis::Z) {
        if (grid.hasMinBoundary(Axis::Z)) {
            for (long j = 0; j < grid.Ny; j++) {
                double physical_y = grid.to_y(j, field.getOffset(), field.getOffsetAxis());
                for (long i = 0; i < grid.Nx; i++) {
                    double physical_x = grid.to_x(i, field.getOffset(), field.getOffsetAxis());
                    dzz(i, j, 0) = (4.0 / 3.0 * field(i, j, 1) - 4.0 * field(i, j, 0) 
                                 + (8.0 / 3.0) *  bc(physical_x, physical_y, 0.0, time)) * mul;
                }
            }
        }

        if (grid.hasMaxBoundary(Axis::Z)) {
            for (long j = 0; j < grid.Ny; j++) {
                for (long i = 0; i < grid.Nx; i++) {
                    dzz(i, j, grid.Nz - 1) = 0.0;
                }
            }
        }
    } else {
        if (grid.hasMinBoundary(Axis::Z)) {
            for (long j = 0; j < grid.Ny; j++) {
                for (long i = 0; i < grid.Nx; i++) {
                    dzz(i, j, 0) = 0.0;
                }
            }
        }

        if (grid.hasMaxBoundary(Axis::Z)) {
            for (long j = 0; j < grid.Ny; j++) {
                double physical_y = grid.to_y(j, field.getOffset(), field.getOffsetAxis());
                for (long i = 0; i < grid.Nx; i++) {
                    double physical_x = grid.to_x(i, field.getOffset(), field.getOffsetAxis());
                    dzz(i, j, grid.Nz - 1) = (4.0 / 3.0 *  field(i, j, grid.Nz - 2) 
                                           - 4.0 *  field(i, j, grid.Nz - 1) 
                                           + (8.0 / 3.0) *  bc(physical_x, physical_y, (grid.Nz - 0.5) * grid.dz, time)) * mul;
                }
            }
        }
    }
}