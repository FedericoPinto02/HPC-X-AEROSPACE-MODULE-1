#pragma once
#include <vector>
#include <cmath>

struct gridDimensions {
    size_t Nx, Ny, Nz;
    double dx, dy, dz;
};

struct velocityField{
   
    std::vector<double> u; 
    std::vector<double> v; 
    std::vector<double> w; 

};


struct PhysicalFields{
    
    size_t u_size, p_size;
    gridDimensions field;  // this naming is a bit of a problem
    velocityField velocity; //velocity
      
    //intermediete velocity
    velocityField xi;
    velocityField eta;
    velocityField zeta;

    std::vector<double> p; //pressure

    //intermediate pressure 
    std::vector<double> psi;
    std::vector<double> phi;
    std::vector<double> c_phi;

    // 2nd derivatives, divergence and laplacian
    std::vector<double> div;
    std::vector<double> laplacian;
    std::vector<double> dxx;
    std::vector<double> dyy;
    std::vector<double> dzz;

    // gradient
    velocityField gradient;
    velocityField exactGradient;

    PhysicalFields(const gridDimensions& dims); 
};

// function to get 3D grid indices
void indices(size_t idx, size_t& i, size_t& j, size_t& k, gridDimensions& myGrid);
size_t index(size_t i, size_t j, size_t k, gridDimensions& myGrid);