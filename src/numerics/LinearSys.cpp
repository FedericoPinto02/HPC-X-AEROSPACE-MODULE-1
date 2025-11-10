#include <algorithm>
#include <cmath>
#include <core/Fields.hpp>
#include <numerics/LinearSys.hpp>
#include <numerics/derivatives.hpp>
#include <vector>

LinearSys::LinearSys(int n, BoundaryType boundaryType)
    : boundaryType(boundaryType), matA(n)

{
    unknownX.resize(n, 0.0);
}


void LinearSys::setRhs(const std::vector<double>& newRhs) {
    if (newRhs.size() != matA.getSize()) {
        throw std::runtime_error(
            "Dimension mismatch: newRhs size does not match system size.");
    }
    
    this->rhsC = newRhs;
}

const std::vector<double>& LinearSys::getSolution() const {
    return this->unknownX;
}

void LinearSys::fillSystemPressure(const Field &phi, const Axis direction) {

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
            "Dimension mismatch: rhsIncomplete must be size n.");
    }


    std::shared_ptr<const Grid> grid = phi.getGrid();
    double d = 0;
    switch (direction) {
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
    std::fill(supdiag.begin(), supdiag.end(), -d);
    std::fill(diag.begin(), diag.end(), 1.0 + 2.0 * d);

    // Boundary conditions
    supdiag.front() = - 2.0 * d;
    diag.back() = 1.0 + d;
}

void LinearSys::fillSystemVelocity(
    const SimulationData simData, const VectorField xi,
    const Axis fieldComponent, const Axis derivativeDirection,
    const size_t iStart, const size_t jStart, const size_t kStart) {

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
            "Dimension mismatch: rhsC must be size n.");
    }
    
    // Need to complete RHS, once onto the switch case!

    std::shared_ptr<const Grid> grid = simData.eta.getGrid();
    double dCoef = 0;
    switch (derivativeDirection) {
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

    for (size_t i = 1; i < matA.getSize() - 1; i++) {
        double k = simData.k.valueWithOffset(iStart, jStart, kStart,
                                                derivativeDirection, i);
        double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
        double gamma = simData.dt * simData.nu * 0.5 / beta;
        subdiag[i - 1] = -gamma * dCoef;
        supdiag[i] = -gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;
    }

    switch (boundaryType)
    {
    case BoundaryType::Normal:

    {
        switch (fieldComponent)
        {
        case Axis::X:
            rhsC.front() =
                simData.bcu(simData.currTime, 0.0, jStart * simData.dy, kStart * simData.dz) -
                grid->dx / 2 *
                    ((simData.bcv(simData.currTime, 0.0, (jStart + 0.5) * simData.dy, kStart * simData.dz) -
                      simData.bcv(simData.currTime, 0.0, (jStart - 0.5) * simData.dy, kStart * simData.dz)) /
                         grid->dy +
                     (simData.bcw(simData.currTime, 0.0, jStart * simData.dy, (kStart + 0.5) * simData.dz) -
                      simData.bcw(simData.currTime, 0.0, jStart * simData.dy, (kStart - 0.5) * simData.dz)) /
                         grid->dz);
            break;

        case Axis::Y:
            rhsC.front() =
                simData.bcv(simData.currTime, iStart * simData.dx, 0.0, kStart * simData.dz) -
                grid->dy / 2 *
                    ((simData.bcu(simData.currTime, (iStart + 0.5) * simData.dx, 0.0, kStart * simData.dz) -
                      simData.bcu(simData.currTime, (iStart - 0.5) * simData.dx, 0.0, kStart * simData.dz)) /
                         grid->dx +
                     (simData.bcw(simData.currTime, iStart * simData.dx, 0.0, (kStart + 0.5) * simData.dz) -
                      simData.bcw(simData.currTime, iStart * simData.dx, 0.0, (kStart - 0.5) * simData.dz)) /
                         grid->dz);
            break;

        case Axis::Z:
            rhsC.front() =
                simData.bcw(simData.currTime, iStart * simData.dx, jStart * simData.dy, 0.0) -
                grid->dz / 2 *
                    ((simData.bcu(simData.currTime, (iStart + 0.5) * simData.dx, jStart * simData.dy, 0.0) -
                      simData.bcu(simData.currTime, (iStart - 0.5) * simData.dx, jStart * simData.dy, 0.0)) /
                         grid->dx +
                     (simData.bcv(simData.currTime, iStart * simData.dx, (jStart + 0.5) * simData.dy, 0.0) -
                      simData.bcv(simData.currTime, iStart * simData.dx, (jStart - 0.5) * simData.dy, 0.0)) /
                         grid->dy);
            break;

        default:
            break;
        }

        diag.back() = 1.0; // the other is aready initialized to zero
        diag.front() = 1.0;

        rhsC.back() = simData.uBoundNew(fieldComponent)
                          .valueWithOffset(iStart, jStart, kStart,
                                           derivativeDirection, matA.getSize() - 1);

        break;
    }

    case BoundaryType::Tangent:
    {
        diag.front() = 1.0;
        double k = simData.k.valueWithOffset(
            iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1);
        double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
        double gamma = simData.dt * simData.nu * 0.5 / beta;
        diag.back() = 1 + 3 * gamma * dCoef;
        subdiag.back() = -gamma * dCoef;

        rhsC.front() = simData.uBoundNew(fieldComponent)
                           .valueWithOffset(iStart, jStart, kStart,
                                            derivativeDirection, 0);

        std::function<double (double t, double x, double y, double z)> bc;   
        double xOff = 0.0, yOff = 0.0, zOff = 0.0;                 
        switch (fieldComponent)
        {
        case Axis::X:
            bc = simData.bcu;
            xOff = 0.5;
            break;

        case Axis::Y:
            bc = simData.bcv;
            yOff = 0.5;
            break;

        case Axis::Z:
            bc = simData.bcw;
            zOff = 0.5;
            break;
        default:
            break;
        }

        switch (derivativeDirection)
        {
        case Axis::X:
            rhsC.back() =
            xi(fieldComponent)
                .valueWithOffset(iStart, jStart, kStart,
                                 derivativeDirection, matA.getSize() - 1) +

            -gamma * dCoef * simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 2) +

            3.0 * gamma * dCoef *
                simData.eta(fieldComponent)
                    .valueWithOffset(iStart, jStart, kStart,
                                     derivativeDirection, matA.getSize() - 1) +
            -2.0 * gamma * dCoef *  // uOld
                bc(simData.currTime - simData.dt, (iStart + matA.getSize() - 1)*simData.dx, (jStart + yOff)*simData.dy, (kStart + zOff)*simData.dz) + 

            2.0 * gamma * dCoef *
                 bc(simData.currTime, (iStart + matA.getSize() - 1)*simData.dx, (jStart + yOff)*simData.dy, (kStart + zOff)*simData.dz);
            break;

        case Axis::Y: 
            rhsC.back() =
                xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1) +
                -gamma * dCoef * simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 2) +
                3.0 * gamma * dCoef * simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1) +
                
                
                -2.0 * gamma * dCoef * bc(simData.currTime - simData.dt, (iStart + xOff)*simData.dx, (jStart + matA.getSize() - 1)*simData.dy, (kStart + zOff)*simData.dz) + 
                
                2.0 * gamma * dCoef *
                    bc(simData.currTime, (iStart + xOff)*simData.dx, (jStart + matA.getSize() - 1)*simData.dy, (kStart + zOff)*simData.dz);
            break;

        case Axis::Z: 
            rhsC.back() =
                xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1) +
                -gamma * dCoef * simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 2) +
                3.0 * gamma * dCoef * simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1) +

                -2.0 * gamma * dCoef * bc(simData.currTime - simData.dt, (iStart + xOff)*simData.dx, (jStart + yOff)*simData.dy, (kStart + matA.getSize() - 1)*simData.dz) + 
                
                2.0 * gamma * dCoef *
                    bc(simData.currTime, (iStart + xOff)*simData.dx, (jStart + yOff)*simData.dy, (kStart + matA.getSize() - 1)*simData.dz);
            break;
        default:
            break;
        }
        
        break;
    }
    default:
        break;
    }
}


void LinearSys::ThomaSolver() {
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
    int size = matA.getSize();

    c[0] /= b[0];
    rhsC[0] /= b[0];
    for (unsigned int i = 1; i < size; i++) {
        b[i] = b[i] - a[i - 1] * c[i - 1];
        rhsC[i] = rhsC[i] - a[i - 1] * rhsC[i - 1];
        if (i < size - 1) // needed because of c dimensions (n-1)
            c[i] /= b[i];
        rhsC[i] /= b[i];
    }

    // Backward Sweep
    unknownX[size - 1] = rhsC[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        unknownX[i] = rhsC[i] - c[i] * unknownX[i + 1];
    }
};
