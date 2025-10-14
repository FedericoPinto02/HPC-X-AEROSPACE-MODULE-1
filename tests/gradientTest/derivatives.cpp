#include "derivatives.hpp"


void computeGradient(PhysicalFields& field){
// compute numerical gradient with finite differences
    size_t i,j,k, idx;
    

// loop over 3D grid indices
    for (size_t i = 1; i < (field.field.Nx-1); i++){
        for (size_t j = 1; j <( field.field.Ny-1); j++){
            for (size_t k = 1; k < (field.field.Nz-1); k++){
                idx =index(i,j,k, field.field);

                field.gradient.u[idx] = (field.p[index(i+1,j,k, field.field)] - field.p[index(i-1,j,k, field.field)]) / (2 * field.field.dx);
                field.gradient.v[idx] = (field.p[index(i,j+1,k, field.field)] - field.p[index(i,j-1,k, field.field)]) / (2 * field.field.dy);
                field.gradient.w[idx] = (field.p[index(i,j,k+1, field.field)] - field.p[index(i,j,k-1, field.field)]) / (2 * field.field.dz);
                
                // write log
                std::cout << ">> position "<<i<<" "<<j<<" "<<k << std::endl;
                std::cout << " pressure "<<field.p[index(i,j,k, field.field)] << std::endl;
                std::cout << " d "<< field.field.dx <<"  "<< field.field.dy<<"  " << field.field.dz  << std::endl;
                std::cout << " gradient "<< field.gradient.u[idx] <<"  "<< field.gradient.v[idx]<<"  " << field.gradient.w[idx]  << std::endl;
                std::cout << " exact gradient "<< field.exactGradient.u[idx]<<"  " << field.exactGradient.v[idx]<<"  " << field.exactGradient.w[idx]  << std::endl;

            }
        }
    }

    }