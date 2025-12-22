#pragma once

#include <cstddef>

#include "core/MpiEnv.hpp"
#include <memory>


/// Enum describing the (possibly staggered) position of some points in a grid.
enum class GridStaggering {
    CELL_CENTERED, FACE_CENTERED
};

/// Shared pointer type for Grid struct.
using GridPtr = std::shared_ptr<const struct Grid>;

/// Struct to hold grid dimensions and spacing information along each direction.
struct Grid {
    //==================================================================================================================
    //--- GLOBAL (DOMAIN) GRID INFO ------------------------------------------------------------------------------------
    //==================================================================================================================
    /// Physical lengths of the global grid in each direction.
    double Lx_glob, Ly_glob, Lz_glob;
    /// Number of grid points in the global grid in each direction.
    size_t Nx_glob, Ny_glob, Nz_glob; // as input
    /// Grid spacing in each direction (constant throughout the domain).
    double dx, dy, dz;

    //==================================================================================================================
    //--- LOCAL (CHUNK) GRID INFO --------------------------------------------------------------------------------------
    //==================================================================================================================
    /// Number of grid points in the local grid chunk in each direction (derived from global info and Cartesian topology).
    size_t Nx, Ny, Nz; // derived, owned by this rank
    /// Index-offset in each dimension of this local grid chunk in the global grid.
    long i_start, j_start, k_start;
    /// Number of halo points along each direction.
    size_t n_halo;

    /**
     * @brief Constructor to initialize the grid in an explicitly sequential environment.
     * @param Nx_g the total number of global points in the X direction
     * @param Ny_g the total number of global points in the Y direction
     * @param Nz_g the total number of global points in the Z direction
     * @param dx_ the grid spacing in the X direction
     * @param dy_ the grid spacing in the Y direction
     * @param dz_ the grid spacing in the Z direction
     */
    Grid(size_t Nx_g, size_t Ny_g, size_t Nz_g,
         double dx_, double dy_, double dz_);

    /**
     * @brief Constructor to initialize the grid in a parallel environment.
     * Automatically decomposes the global domain based on MpiEnv topology.
     * @param Nx_g the total number of global points in the X direction
     * @param Ny_g the total number of global points in the Y direction
     * @param Nz_g the total number of global points in the Z direction
     * @param dx_ the grid spacing in the X direction
     * @param dy_ the grid spacing in the Y direction
     * @param dz_ the grid spacing in the Z direction
     * @param env the MPI Environment containing topology information
     * @param n_halo the number of halo points along each direction (default is 1)
     */
    Grid(size_t Nx_g, size_t Ny_g, size_t Nz_g,
         double dx_, double dy_, double dz_,
         const MpiEnv &env,
         size_t n_halo = 1);

    /**
     * @brief Decompose the grid along a target dimension, identifying local information of the chunk owned by a process
     * with a specified coordinate in a Cartesian MPI topology (parallel environment).
     * The domain decomposition is done in a balanced manner, distributing any remainder points to the first few.
     * @param N_glob the total number of global points in the target dimension
     * @param n_procs the number of processes along the target dimension
     * @param coord the coordinate of the current process along the target dimension
     * @return a pair containing the number of local points for the target dimension
     *          and the starting index (offset) for that dimension in the global grid
     */
    static std::pair<size_t, size_t> decompose1D(size_t N_glob, int n_procs, int proc_coord);

    /**
     * @brief Getter for the total number of global grid points.
     * @return the total number of global grid points
     */
    [[nodiscard]] inline constexpr size_t globalSize() const { return Nx_glob * Ny_glob * Nz_glob; }

    /**
     * @brief Getter for the total number of local grid points.
     * @return the total number of local grid points
     */
    [[nodiscard]] inline constexpr size_t size() const { return Nx * Ny * Nz; }

    /**
     * @brief Converts a grid index to a physical coordinate along the X-axis, considering staggering.
     * @param i the x-index
     * @param offset the staggering offset
     * @param offsetAxis the axis where to apply the offset
     * @return the physical x-coordinate
     */
    [[nodiscard]] inline double to_x(long i, GridStaggering offset, Axis offsetAxis) const {
        const auto i_glob = (double) (i_start + i);
        const auto i_offset = (offset == GridStaggering::FACE_CENTERED && offsetAxis == Axis::X) ? 0.5 : 0.0;
        return (i_glob + i_offset) * dx;
    }

    /**
     * @brief Converts a grid index to a physical coordinate along the Y-axis, considering staggering.
     * @param j the y-index
     * @param offset the staggering offset
     * @param offsetAxis the axis where to apply the offset
     * @return the physical y-coordinate
     */
    [[nodiscard]] inline double to_y(long j, GridStaggering offset, Axis offsetAxis) const {
        const auto j_glob = (double) (j_start + j);
        const auto j_offset = (offset == GridStaggering::FACE_CENTERED && offsetAxis == Axis::Y) ? 0.5 : 0.0;
        return (j_glob + j_offset) * dy;
    }

    /**
     * @brief Converts a grid index to a physical coordinate along the Z-axis, considering staggering.
     * @param k the z-index
     * @param offset the staggering offset
     * @param offsetAxis the axis where to apply the offset
     * @return the physical z-coordinate
     */
    [[nodiscard]] inline double to_z(long k, GridStaggering offset, Axis offsetAxis) const {
        const auto k_glob = (double) (k_start + k);
        const auto k_offset = (offset == GridStaggering::FACE_CENTERED && offsetAxis == Axis::Z) ? 0.5 : 0.0;
        return (k_glob + k_offset) * dz;
    }

    /**
     * @brief Checks if the local grid chunk has the minimum boundary at the specified axis.
     * @param axis the axis to check (Axis::X, Axis::Y, or Axis::Z)
     * @return true if the local grid chunk has the minimum boundary at the specified axis, false otherwise
     */
    [[nodiscard]] inline bool hasMinBoundary(Axis axis) const {
        if (axis == Axis::X) return i_start == 0;
        if (axis == Axis::Y) return j_start == 0;
        return k_start == 0;
    }

    /**
     * @brief Checks if the local grid chunk has the maximum boundary at the specified axis.
     * @param axis the axis to check (Axis::X, Axis::Y, or Axis::Z)
     * @return true if the local grid chunk has the maximum boundary at the specified axis, false otherwise
     */
    [[nodiscard]] inline bool hasMaxBoundary(Axis axis) const {
        if (axis == Axis::X) return (i_start + Nx) == Nx_glob;
        if (axis == Axis::Y) return (j_start + Ny) == Ny_glob;
        return (k_start + Nz) == Nz_glob;
    }
};
