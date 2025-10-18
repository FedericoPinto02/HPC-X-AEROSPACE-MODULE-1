#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"
#include <stdexcept>

void Derivatives::computeDx(const Field &field, Field &dx) {

    // safety checks
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dx.getGrid()) {
        throw std::runtime_error("output field (dx) has null grid.");
    }
    if (field.getGrid() != dx.getGrid()) {
        throw std::runtime_error("field and dx must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 0.5/grid->dx;

    if (Nx < 3) {
        throw std::runtime_error("grid size Nx must be at least 3 to compute central derivative.");
    }
    if (grid->dx <= 0.0) {
        throw std::runtime_error("grid spacing dx must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=0 ; k<Nz ; k++)
        for (size_t j=0 ; j<Ny ; j++)
            for (size_t i=1 ; i<(Nx-1) ; i++)
                {
                    double dfdx = 0.0;
                    dfdx = (field(i + 1, j, k) - field(i - 1, j, k))*mul;
                    dx(i,j,k) = dfdx;
                }
}

void Derivatives::computeDy(const Field &field, Field &dy) {

    // safety checks
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dy.getGrid()) {
        throw std::runtime_error("output field (dy) has null grid.");
    }
    if (field.getGrid() != dy.getGrid()) {
        throw std::runtime_error("field and dx must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 0.5/grid->dy;

    if (Ny < 3) {
        throw std::runtime_error("grid size Ny must be at least 3 to compute central derivative.");
    }
    if (grid->dy <= 0.0) {
        throw std::runtime_error("grid spacing dx must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=0 ; k<Nz ; k++)
        for (size_t i=0 ; i<Nx ; i++)
            for (size_t j=1 ; j<(Ny-1) ; j++)
                {
                    double dfdy = 0.0;
                    dfdy = (field(i, j + 1, k) - field(i, j - 1, k)) * mul;
                    dy(i,j,k) = dfdy;
                }
}