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
    case Axis::X:
        d = 1 / grid->dx / grid->dx;
        break;
    case Axis::Y:
        d = 1 / grid->dy / grid->dy;
        break;
    case Axis::Z:
        d = 1 / grid->dz / grid->dz;
        break;
    default:
        break;
    }
    std::fill(subdiag.begin(), subdiag.end(), -d);
    subdiag.back() *= 2;
    std::fill(supdiag.begin(), supdiag.end(), -d);
    supdiag.front() *= 2;
    std::fill(diag.begin(), diag.end(), 1+2*d);
}


 void LinearSys::fillSystemVelocity(Field& porosity, std::vector<double>& rhsIncomplete, 
        VectorField& eta, VectorField& xi, VectorField& uBoundNew, VectorField& uBoundOld,
        const Axis fieldComponent, const Axis derivativeDirection, const size_t iStart, const size_t jStart, const size_t kStart)
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

    std::shared_ptr<const Grid> grid = eta.getGrid();
    double dCoef = 0;
    switch (derivativeDirection)
    {
    case Axis::X:
        dCoef = 1 / grid->dx / grid->dx;
        break;
    case Axis::Y:
        dCoef = 1 / grid->dy / grid->dy;
        break;
    case Axis::Z:
        dCoef = 1 / grid->dz / grid->dz;
        break;
    default:
        break;
    }

    for (size_t i=1 ; i<matA.getSize()-1 ; i++)
    {
        double gamma = porosity.valueWithOffset(iStart, jStart, kStart, derivativeDirection, i);
        subdiag[i-1] = - gamma * dCoef;
        supdiag[i] = - gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;
    }
    
    switch (boundaryType)
    {
    case BoundaryType::Normal:
        
        switch (fieldComponent)
        {
        case Axis::X:  // it's a dirichlet on u component
            rhsC.front() = uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart, Axis::X,  0) 
                            - grid->dx /2 * ( (uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart, Axis::Y,  1)
                                            - uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart, Axis::Y,  -1)) / grid->dy +
                                            (uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart, Axis::Z,  1)
                                            - uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart, Axis::Z,  -1)) / grid->dz
                            );
            break;
        case Axis::Y:
            rhsC.front() = uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart, Axis::Y,  0) 
                            - grid->dy /2 * ( (uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart, Axis::X,  1)
                                            - uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart, Axis::X,  -1)) / grid->dx +
                                            (uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart, Axis::Z,  1)
                                            - uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart, Axis::Z,  -1)) / grid->dz
                            );
            break;
        case Axis::Z:
            rhsC.front() = uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart, Axis::Z,  0) 
                            - grid->dz /2 * ( (uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart, Axis::Y,  1)
                                                - uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart, Axis::Y,  -1)) / grid->dy +
                                                (uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart, Axis::X,  1)
                                                - uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart, Axis::X,  -1)) / grid->dx
                            );
            break;
        default:
            break;
        }
        
        diag.back() = 1.0;  // the other is aready initialized to zero
        

        rhsC.back() = uBoundNew(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize()); 
        break;

    
    case BoundaryType::Tangent:
        diag.front() = 1.0;
        double gamma = porosity.valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize());
        diag.back() = 1 + 3 * gamma * dCoef;
        subdiag.back() = - gamma * dCoef;

        rhsC.front() = uBoundNew(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,0); 

        rhsC.back() = 2*gamma*dCoef *uBoundNew(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,matA.getSize()) + 
                        2*gamma*dCoef *uBoundOld(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,matA.getSize()) + 
                        gamma * dCoef * eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,matA.getSize()-1) +
                        -3* gamma * dCoef * eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,matA.getSize()) +
                         gamma * dCoef * xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,matA.getSize());
        break;
    default:
        break;
    }
        
}

void LinearSys::ThomaSolver(){
    /*
    Solves the tridiagonal system Ax = f using the Thomas algorithm.
    a: sub-diagonal (a[0] is not used)
    b: main diagonal
    c: super-diagonal (c[n-1] is not used)
    f: right-hand side
    x: solution vector (to be computed)
    All vectors are of size n, except a and c which are of size n-1.
    A compatibility check is performed to ensure the dimensions are correct.
    */

    // Compatibility check

    std::vector<double> b = matA.getDiag(0);
    std::vector<double> a = matA.getDiag(-1);
    std::vector<double> c = matA.getDiag(1);
   
    c[0] /= b[0];
    rhsC[0] /= b[0];
    for ( unsigned int i = 1; i < matA.getSize() ; i++)
    {
        b[i] = b[i] - a[i-1] * c[i-1];
        rhsC[i] = rhsC[i] - a[i-1] * rhsC[i-1];
        if (i < matA.getSize() - 1)                  // needed because of c dimensions (n-1)
            c[i] /= b[i];
        rhsC[i] /= b[i];
    }

    // Backward Sweep 
    unknownX[matA.getSize()-1] = rhsC[matA.getSize()-1];
    for ( int i = matA.getSize() - 2; i >= 0; i-- )
    {
        unknownX[i] = rhsC[i] - c[i] * unknownX[i+1];
    }
};



