#include "physicalField.hpp"



size_t index(size_t i, size_t j, size_t k, gridDimensions& myGrid)  {
    return i + myGrid.Nx * (j + myGrid.Ny * k);
}

void indices(size_t idx, size_t& i, size_t& j, size_t& k, gridDimensions& myGrid) {
    k = idx / (myGrid.Nx * myGrid.Ny);
    j = (idx / myGrid.Nx) % myGrid.Ny;
    i = idx % myGrid.Nx;
}

PhysicalFields::PhysicalFields(const gridDimensions& dims){

    field = dims;
    u_size = field.Nx * field.Ny * field.Nz;
    p_size = field.Nx * field.Ny * field.Nz;

    velocity.u.resize(u_size, 0.0);
    velocity.v.resize(u_size, 0.0);
    velocity.w.resize(u_size, 0.0);

    gradient.u.resize(u_size, 0.0);
    gradient.v.resize(u_size, 0.0);
    gradient.w.resize(u_size, 0.0);

    exactGradient.u.resize(u_size, 0.0);
    exactGradient.v.resize(u_size, 0.0);
    exactGradient.w.resize(u_size, 0.0);


    p.resize(p_size, 0.0);

    div.resize(p_size, 0.0);


}