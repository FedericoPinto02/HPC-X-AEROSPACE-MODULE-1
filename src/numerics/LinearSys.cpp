 #include <vector>
#include <numerics/LinearSys.hpp>
#include <cmath>
#include <numerics/derivatives.hpp>
#include <algorithm>
#include <core/Fields.hpp>


LinearSys::LinearSys(int n, BoundaryType boundaryType) 
    : 
     boundaryType (boundaryType),
     matA (n)

{
    rhsC.resize(n, 0.0);
    unknownX.resize(n, 0.0);
}

void LinearSys::fillSystemPressure(std::vector<double>& rhsIncomplete, Field& phi, const Axis direction)
{

    std::vector<double>& diag = matA.getDiag(0);
    std::vector<double>& subdiag = matA.getDiag(-1);
    std::vector<double>& supdiag = matA.getDiag(1);

    if (rhsIncomplete.size() != matA.getSize() - 2) {
        throw std::runtime_error("Dimension mismatch: rhsIncomplete must be size - 2.");
    }
    // Assign incomplete RHS 
    std::copy(
        rhsIncomplete.begin(), 
        rhsIncomplete.end(),
        rhsC.begin() + 1       
    );
     // No need to complete RHS, already set to 0.0
    std::shared_ptr<const Grid> grid = phi.getGrid();
    double d = 0;
    switch (direction)
    {
    case Axis::x:
        d = 1 / grid->dx / grid->dx;
        break;
    case Axis::y:
        d = 1 / grid->dy / grid->dy;
        break;
    case Axis::z:
        d = 1 / grid->dz / grid->dz;
        break;
    default:
        break;
    }
    std::fill(subdiag.begin(), subdiag.end(), -d);
    subdiag.front() *= 2;
    subdiag.back() *= 2;
    std::fill(supdiag.begin(), supdiag.end(), -d);
    supdiag.front() *= 2;
    supdiag.back() *= 2;
    std::fill(diag.begin(), diag.end(), 1+2*d);
}

// fix iStart please
 void LinearSys::fillSystemVelocity(Field& porosity, std::vector<double>& rhsIncomplete, 
        VectorField& etaNew, VectorField& eta, VectorField& xi,
        const Axis direction, const size_t iStart, const size_t jStart, const size_t kStart)
{

    std::vector<double>& diag = matA.getDiag(0);
    std::vector<double>& subdiag = matA.getDiag(-1);
    std::vector<double>& supdiag = matA.getDiag(1);

    if (rhsIncomplete.size() != matA.getSize() - 2) {
        throw std::runtime_error("Dimension mismatch: rhsIncomplete must be size - 2.");
    }
    // Assign incomplete RHS 
    std::copy(
        rhsIncomplete.begin(), 
        rhsIncomplete.end(),
        rhsC.begin() + 1       
    );
    // Need to complete RHS, once onto the switch case!

    std::shared_ptr<const Grid> grid = etaNew.getGrid();
    double d = 0;
    switch (direction)
    {
    case Axis::x:
        d = 1 / grid->dx / grid->dx;
        break;
    case Axis::y:
        d = 1 / grid->dy / grid->dy;
        break;
    case Axis::z:
        d = 1 / grid->dz / grid->dz;
        break;
    default:
        break;
    }

    for (size_t i=1 ; i<matA.getSize() ; i++)
    {
        double gamma = porosity(iStart, jStart, kStart, direction, i);
        subdiag[i-2] = - gamma * d;
        supdiag[i] = - gamma * d;
        diag[i] = 1 + 2 * gamma * d;
    }
    
    switch (boundaryType)
    {
        // TO DO: write b.c. for RHS
    case BoundaryType::Normal:
        /** 
        switch (direction)
        {
        case Axis::x:
            rhsC = 1 / grid->dx / grid->dx;
            break;
        case Axis::y:
            d = 1 / grid->dy / grid->dy;
            break;
        case Axis::z:
            d = 1 / grid->dz / grid->dz;
            break;
        default:
            break;
        }
        */
        diag.back() = 1.0;  // the other is aready initialized to zero
        rhsC.back() = eta
        break;
    
    case BoundaryType::Tangent:
        diag.front() = 1.0;
        double gamma = porosity(iStart, jStart, kStart, direction, matA.getSize());
        diag.back() = 1 + 3 * gamma * d;
        subdiag.back() = - gamma * d;

    default:
        break;
    }
        
    

}



