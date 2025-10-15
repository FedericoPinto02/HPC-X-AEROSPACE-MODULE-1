#include "Derivatives.hpp"

Derivatives::Derivatives(double dx_, double dy_, double dz_)
    : dx(dx_), dy(dy_), dz(dz_){}

std::vector<double> Derivatives::derive( std::vector<double>& v, double ds ) const
{
    size_t n = v.size();
    std::vector<double> dv(n,0.0);
    double step = 0.5/ds;

    for (size_t i = 1; i<n-1 ;i++)
        dv[i] = (v[i+1] - v[i-1])*step;

    return dv;
}