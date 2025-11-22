#ifndef NSBSOLVER_MESH_HPP
#define NSBSOLVER_MESH_HPP

#include <cstddef>

/**
 * @brief Enum class representing the three coordinate axes.
 */
enum class Axis {
    X = 0, Y = 1, Z = 2
};

const size_t AXIS_COUNT = 3;

/**
 * @brief Enum describing the (possibly staggered) position of some points in a grid.
 */
enum class GridStaggering {
    CELL_CENTERED, FACE_CENTERED
};

/**
* @brief Struct to hold grid dimensions and spacing information along each direction.
*/
struct Grid {
    /// Physical size of the grid in each direction.
    double Lx, Ly, Lz;
    /// Number of grid points in each direction.
    size_t Nx, Ny, Nz;
    /// Grid spacing in each direction (derived from physical sizes and number of grid points).
    double dx, dy, dz;

    Grid() {
        Nx = Ny = Nz = 10;
        dx = dy = dz = 1.0;
        Lx = dx * (static_cast<double>(Nx) + 0.5);
        Ly = dy * (static_cast<double>(Ny) + 0.5);
        Lz = dz * (static_cast<double>(Nz) + 0.5);
    }

    /**
     * @brief Constructor to initialize the grid with given physical sizes and number of points.
     * @param Nx_ the number of grid points in the X direction
     * @param Ny_ the number of grid points in the Y direction
     * @param Nz_ the number of grid points in the Z direction
     * @param dx_ the grid spacing in the X direction
     * @param dy_ the grid spacing in the Y direction
     * @param dz_ the grid spacing in the Z direction
     */
    Grid(size_t Nx_, size_t Ny_, size_t Nz_, double dx_, double dy_, double dz_)
            : Nx(Nx_), Ny(Ny_), Nz(Nz_), dx(dx_), dy(dy_), dz(dz_) {
        if (Nx <= 0 || Ny <= 0 || Nz <= 0) {
            throw std::runtime_error("Grid dimensions must be positive.");
        }
        if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
            throw std::runtime_error("Grid spacing must be positive.");
        }
        Lx = dx * (static_cast<double>(Nx) + 0.5);
        Ly = dy * (static_cast<double>(Ny) + 0.5);
        Lz = dz * (static_cast<double>(Nz) + 0.5);
    }


    /**
     * @brief Getter for the total number of grid points.
     * @return the total number of grid points
     */
    [[nodiscard]] constexpr size_t size() { return Nx * Ny * Nz; }

    /**
     * @overload
     */
    [[nodiscard]] constexpr size_t size() const { return Nx * Ny * Nz; }

    /**
     * @brief Converts a grid index to a physical coordinate along the X-axis, considering staggering.
     * @param i the x-index
     * @param offset the staggering offset
     * @param offsetAxis the axis where to apply the offset
     * @return the physical x-coordinate
     */
    [[nodiscard]] inline double to_x(double i, GridStaggering offset, Axis offsetAxis) const {
        const auto xOffset = (offset == GridStaggering::FACE_CENTERED && offsetAxis == Axis::X) ? 0.5 : 0.0;
        return (i + xOffset) * dx;
    }

    /**
     * @brief Converts a grid index to a physical coordinate along the Y-axis, considering staggering.
     * @param j the y-index
     * @param offset the staggering offset
     * @param offsetAxis the axis where to apply the offset
     * @return the physical y-coordinate
     */
    [[nodiscard]] inline double to_y(double j, GridStaggering offset, Axis offsetAxis) const {
        const auto yOffset = (offset == GridStaggering::FACE_CENTERED && offsetAxis == Axis::Y) ? 0.5 : 0.0;
        return (j + yOffset) * dy;
    }

    /**
     * @brief Converts a grid index to a physical coordinate along the Z-axis, considering staggering.
     * @param k the z-index
     * @param offset the staggering offset
     * @param offsetAxis the axis where to apply the offset
     * @return the physical z-coordinate
     */
    [[nodiscard]] inline double to_z(double k, GridStaggering offset, Axis offsetAxis) const {
        const auto zOffset = (offset == GridStaggering::FACE_CENTERED && offsetAxis == Axis::Z) ? 0.5 : 0.0;
        return (k + zOffset) * dz;
    }
};

#endif // NSBSOLVER_MESH_HPP
