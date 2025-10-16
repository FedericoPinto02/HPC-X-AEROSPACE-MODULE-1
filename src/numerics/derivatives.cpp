#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"

void Derivatives::computeDx(const Field &field, Field &dx) {
std::shared_ptr<const Grid> grid = field.getGrid();
const size_t Nx = grid->Nx;
const size_t Ny = grid->Ny;
const size_t Nz = grid->Nz;
const double mul = 0.5/grid->dx;

for (size_t k=0 ; k<Nz ; k++)
    for (size_t j=0 ; j<Ny ; j++)
        for (size_t i=1 ; i<(Nx-1) ; i++)
            {
                double dfdx = 0.0;
                dfdx = (field(i + 1, j, k) - field(i - 1, j, k))*mul;
                dx(i,j,k) = dfdx;
            }
}