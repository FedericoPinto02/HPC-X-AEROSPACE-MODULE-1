#include "core/Grid.hpp"

Grid::Grid(size_t Nx_g, size_t Ny_g, size_t Nz_g,
           double dx_, double dy_, double dz_)
        : Nx_glob(Nx_g), Ny_glob(Ny_g), Nz_glob(Nz_g),
          Nx(Nx_g), Ny(Ny_g), Nz(Nz_g),
          dx(dx_), dy(dy_), dz(dz_),
          i_start(0L), j_start(0L), k_start(0L),
          n_halo(1) {
    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        throw std::runtime_error("Grid spacing must be positive.");
    }

    // Compute global physical dimensions
    Lx_glob = dx * (static_cast<double>(Nx_glob) - 0.5);
    Ly_glob = dy * (static_cast<double>(Ny_glob) - 0.5);
    Lz_glob = dz * (static_cast<double>(Nz_glob) - 0.5);
}

Grid::Grid(size_t Nx_g, size_t Ny_g, size_t Nz_g,
           double dx_, double dy_, double dz_,
           const MpiEnv &env,
           size_t n_halo_)
        : Nx_glob(Nx_g), Ny_glob(Ny_g), Nz_glob(Nz_g),
          dx(dx_), dy(dy_), dz(dz_),
          n_halo(n_halo_) {
    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        throw std::runtime_error("Grid spacing must be positive.");
    }

    // Compute global physical dimensions
    Lx_glob = dx * (static_cast<double>(Nx_glob) - 0.5);
    Ly_glob = dy * (static_cast<double>(Ny_glob) - 0.5);
    Lz_glob = dz * (static_cast<double>(Nz_glob) - 0.5);

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
    // Safety check for trivial cases
    if (n_procs <= 1) return {N_glob, 0};
    if (N_glob < static_cast<size_t>(n_procs) + 1) {
        throw std::runtime_error("Grid too small for the requested number of processors (Overlapping Decomposition).");
    }

    // --- OVERLAPPING STRATEGY ---
    // 1. Decompose the INTERVALS (edges) between points.
    //    If there are N points, there are N-1 intervals.
    size_t n_intervals = N_glob - 1;

    size_t base = n_intervals / n_procs;
    size_t remainder = n_intervals % n_procs;

    size_t my_intervals = base;
    size_t start_interval = proc_coord * base;

    // Distribute remainder intervals
    if (static_cast<size_t>(proc_coord) < remainder) {
        my_intervals++;
        start_interval += proc_coord;
    } else {
        start_interval += remainder;
    }

    // 2. Convert back to POINTS.
    //    A segment of M intervals contains M+1 points (vertices).
    //    The start index matches the start interval index.
    size_t proc_count = my_intervals + 1;
    size_t proc_start = start_interval;

    return {proc_count, proc_start};
}
