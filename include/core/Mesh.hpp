#ifndef NSBSOLVER_MESH_HPP
#define NSBSOLVER_MESH_HPP

#include <cstddef>


 * @brief Struct to hold grid dimensions and spacing information along each direction.
 */
struct Grid {
    size_t Nx, Ny, Nz; // Number of grid points in each direction
    double dx, dy, dz; // Grid spacing in each direction

    Grid(size_t Nx, size_t Ny, size_t Nz, double dx, double dy, double dz)
            : Nx(Nx), Ny(Ny), Nz(Nz), dx(dx), dy(dy), dz(dz) {};

    /**
     * @brief Getter for the total number of grid points.
     * @return the total number of grid points
     */
    size_t size() { return Nx * Ny * Nz; }

    /**
     * @overload
     */
    size_t size() const { return Nx * Ny * Nz; }
};
// todo - TBD: Mesh class ?? methods for grid generation ??

#endif // NSBSOLVER_MESH_HPP
