#include "core/Grid.hpp"

Grid::Grid(size_t Nx_g, size_t Ny_g, size_t Nz_g,
           double dx_, double dy_, double dz_)
        : Nx_glob(Nx_g), Ny_glob(Ny_g), Nz_glob(Nz_g),
          Nx(Nx_g), Ny(Ny_g), Nz(Nz_g),
          dx(dx_), dy(dy_), dz(dz_),
          i_start(0), j_start(0), k_start(0) {
    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        throw std::runtime_error("Grid spacing must be positive.");
    }

    // Compute global physical dimensions
    Lx_glob = dx * (static_cast<double>(Nx_glob) + 0.5);
    Ly_glob = dy * (static_cast<double>(Ny_glob) + 0.5);
    Lz_glob = dz * (static_cast<double>(Nz_glob) + 0.5);
}

Grid::Grid(size_t Nx_g, size_t Ny_g, size_t Nz_g,
           double dx_, double dy_, double dz_,
           const MpiEnv &env)
        : Nx_glob(Nx_g), Ny_glob(Ny_g), Nz_glob(Nz_g),
          dx(dx_), dy(dy_), dz(dz_) {
    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        throw std::runtime_error("Grid spacing must be positive.");
    }

    // Compute global physical dimensions
    Lx_glob = dx * (static_cast<double>(Nx_glob) + 0.5);
    Ly_glob = dy * (static_cast<double>(Ny_glob) + 0.5);
    Lz_glob = dz * (static_cast<double>(Nz_glob) + 0.5);

    // Compute own chunk information based on MPI topology
    auto dims = env.dims();     // e.g., {2, 2, 1}
    auto coords = env.coords(); // e.g., {0, 1, 0}
    std::tie(Nx, i_start) = decompose1D(Nx_glob, dims[0], coords[0]);
    std::tie(Ny, j_start) = decompose1D(Ny_glob, dims[1], coords[1]);
    std::tie(Nz, k_start) = decompose1D(Nz_glob, dims[2], coords[2]);
}

std::pair<size_t, size_t> Grid::decompose1D(
        size_t N_glob,
        int n_procs,
        int proc_coord
) {
    if (n_procs == 0) return {N_glob, 0}; // Safety for non-MPI runs

    size_t base = N_glob / n_procs;
    size_t remainder = N_glob % n_procs;

    size_t proc_count = base;
    size_t proc_start = proc_coord * base;

    // Distribute the remainder points to the first few processors
    if (static_cast<size_t>(proc_coord) < remainder) {
        proc_count++;
        proc_start += proc_coord; // I have 'coord' extra points before me
    } else {
        proc_start += remainder; // All remainders are before me
    }

    return {proc_count, proc_start};
}

bool Grid::hasMinBoundary(Axis axis) const {
    if (axis == Axis::X) return i_start == 0;
    if (axis == Axis::Y) return j_start == 0;
    return k_start == 0;
}

bool Grid::hasMaxBoundary(Axis axis) const {
    if (axis == Axis::X) return (i_start + Nx) == Nx_glob;
    if (axis == Axis::Y) return (j_start + Ny) == Ny_glob;
    return (k_start + Nz) == Nz_glob;
}
