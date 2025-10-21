#ifndef NSBSOLVER_MESH_HPP
#define NSBSOLVER_MESH_HPP

#include <cstddef>


/**
 * @brief Struct to hold grid dimensions and spacing information along each direction.
 */
struct Grid {
    size_t Nx, Ny, Nz;
    double dx, dy, dz;

    Grid(size_t Nx, size_t Ny, size_t Nz, double dx, double dy, double dz)
        : Nx(Nx), Ny(Ny), Nz(Nz), dx(dx), dy(dy), dz(dz) {}

    size_t size() const { return Nx * Ny * Nz; }
};

// todo - TBD: Mesh class ?? methods for grid generation ??

#endif // NSBSOLVER_MESH_HPP