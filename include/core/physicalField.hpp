#include <vector>
#include <cmath>

struct gridDimensions {
    size_t Nx, Ny, Nz;
};

struct velocityField{
   
    std::vector<double> u; 
    std::vector<double> v; 
    std::vector<double> w; 

};


struct PhysicalFields{
    
    size_t u_size, p_size;
    gridDimensions field;
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
};

