#include "Derivatives.hpp"

Derivatives::Derivatives(double dx_, double dy_, double dz_, double nx_, double ny_, double nz_)  
: dx(dx_), dy(dy_), dz(dz_), nx(nx_), ny(ny_), nz(nz_){}

std::vector<double> Derivatives::derive_x( std::vector<double>& v) const
{
    std::vector<double> dv(v.size(),0.0);
    double step = 0.5/dx;
    size_t index = 0;

    for (size_t j = 1; j<ny-1 ;j++)
        for (size_t k = 1; k<nz-1; k++)
            for (size_t i = 1; i<nx-1; i++) {
                index = i + nx*(j + ny*k);
                dv[index] = (v[index+1] - v[index-1])*step;
            }

    return dv;
}