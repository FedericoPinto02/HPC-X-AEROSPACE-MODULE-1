#include <iostream>
#include <cmath>
#include "physicalField.hpp"
#include "derivatives.hpp"




int main(){

    gridDimensions myGrid = {10, 10, 10, 0.001, 0.001, 0.001}; // Nx, ... dx,...
    PhysicalFields myField(myGrid);
    size_t i,j,k; 
    double pi = 3.141; // shitty.pi
    double w1 = 0.0123, w2 = 0.4423, w3 = 0.7801; // known signal periods

    for (size_t idx = 0; idx < myField.p_size; idx++){
        indices(idx, i, j, k, myGrid);
        // pressure
        myField.p[idx] = std::sin(i * pi * w1 / myGrid.Nx) + std::sin(j * pi * w2/ myGrid.Ny) + std::sin(k * pi * w3/myGrid.Nz);
        // exact gradient
        myField.exactGradient.u[idx] = (pi*w1)/(myGrid.Nx* myGrid.dx) * std::cos(i * pi * w1/ myGrid.Nx);
        myField.exactGradient.v[idx] = (pi*w2)/(myGrid.Ny* myGrid.dy) * std::cos(j * pi * w2/ myGrid.Ny);
        myField.exactGradient.w[idx] = (pi*w3)/(myGrid.Nz* myGrid.dz) * std::cos(k * pi * w3 /myGrid.Nz);
    }

    computeGradient(myField);

    // how to derive 
    // sin(x/ X *pi*w1)  ----> 1/X *pi*w1 *cos(x/X *pi *w1)
    // sin((i*dx )/ (Nx * dx) *pi*w1)  ----> 1/X *pi*w1 *cos(x/X *pi *w1)
 



 

    return 0;
}
