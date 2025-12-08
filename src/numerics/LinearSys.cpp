#include <algorithm>
#include <cmath>
#include <core/Fields.hpp>
#include <numerics/LinearSys.hpp>
#include <numerics/derivatives.hpp>
#include <vector>

LinearSys::LinearSys(int n)
        : matA(n) {
    unknownX.resize(n, 0.0);
}


void LinearSys::setRhs(const std::vector<double> &newRhs) {
    if (newRhs.size() != matA.getSize()) {
        throw std::runtime_error(
                "Dimension mismatch: newRhs size does not match system size.");
    }

    this->rhsC = newRhs;
}

const std::vector<double> &LinearSys::getSolution() const {
    return this->unknownX;
}

void LinearSys::fillSystemPressure(const GridPtr &grid, const Axis direction) {

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
                "Dimension mismatch: rhsIncomplete must be size n.");
    }

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
    supdiag.front() = -2.0 * d;
    diag.back() = 1.0 + d;
}

void LinearSys::fillSystemVelocity(
        const SimulationData &simData, const VectorField &eta, const VectorField &xi,
        const Axis fieldComponent, const Axis derivativeDirection,
        const size_t iStart, const size_t jStart, const size_t kStart) {

    boundaryType = (fieldComponent == derivativeDirection)
               ? BoundaryType::Normal
               : BoundaryType::Tangent;

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
                "Dimension mismatch: rhsC must be size n.");
    }

    // Need to complete RHS, once onto the switch case!

    const Grid &grid = eta.getGrid();
    double dCoef = 0;
    switch (derivativeDirection) {
        case Axis::X:
            dCoef = 1 / grid.dx / grid.dx;
            break;
        case Axis::Y:
            dCoef = 1 / grid.dy / grid.dy;
            break;
        case Axis::Z:
            dCoef = 1 / grid.dz / grid.dz;
            break;
        default:
            break;
    }

    for (size_t i = 1; i < matA.getSize() - 1; i++) {
        double k = simData.k(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                             derivativeDirection, i);
        double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
        double gamma = simData.dt * simData.nu * 0.5 / beta;
        subdiag[i - 1] = -gamma * dCoef;
        supdiag[i] = -gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;
    }

    switch (boundaryType) {
        case BoundaryType::Normal: {
            double k = simData.k(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, 0);
            double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
            double gamma = simData.dt * simData.nu * 0.5 / beta;
            diag.front() = 1 + 4.0 * gamma * dCoef;
            supdiag.front() = -4.0 / 3.0  * gamma * dCoef;
            diag.back() = 1.0; // the other is aready initialized to zero

            std::function<double(double t, double x, double y, double z)> bc;
            double xOff = 0.0, yOff = 0.0, zOff = 0.0;
            switch (fieldComponent) {  // if normal fieldComponent == derivativeDirection
                case Axis::X:
                    bc = simData.bcu;
                    xOff = 0.5;
                    rhsC.front() =
                            xi(fieldComponent)(0, jStart, kStart) 
                            + 4.0 /3.0 * gamma * dCoef *  eta(fieldComponent)(1,jStart,kStart)
                            - 4.0 * gamma * dCoef * eta(fieldComponent)(0, jStart, kStart) 
                            + 8.0/3.0 * gamma * dCoef *  (
                            bc(0.0, jStart * grid.dy, kStart * grid.dz, simData.currTime - simData.dt) +
                            bc(0.0, jStart * grid.dy, kStart * grid.dz, simData.currTime) );

                    rhsC.back() = bc((iStart + 0.5 + matA.getSize() - 1) * grid.dx, jStart * grid.dy, kStart * grid.dz, simData.currTime);
                    break;

                case Axis::Y:
                    bc = simData.bcv;
                    yOff = 0.5;
                    rhsC.front() =
                            xi(fieldComponent)(iStart, 0, kStart) 
                            + 4.0 /3.0 * gamma * dCoef *  eta(fieldComponent)(iStart, 1, kStart)
                            - 4.0 * gamma * dCoef * eta(fieldComponent)(iStart, 0, kStart) 
                            + 8.0/3.0 * gamma * dCoef *  (
                            bc(iStart * grid.dx, 0.0, kStart * grid.dz, simData.currTime - simData.dt) +
                            bc(iStart * grid.dx, 0.0, kStart * grid.dz, simData.currTime) );

                    rhsC.back() = bc(iStart * grid.dx, (jStart + 0.5 + matA.getSize() - 1) * grid.dy, kStart * grid.dz, simData.currTime);
                    break;

                case Axis::Z:
                    bc = simData.bcw;
                    zOff = 0.5;
                    rhsC.front() =
                            xi(fieldComponent)(iStart, jStart, 0) 
                            + 4.0 /3.0 * gamma * dCoef *  eta(fieldComponent)(iStart,jStart,1)
                            - 4.0 * gamma * dCoef * eta(fieldComponent)(iStart, jStart, 0) 
                            + 8.0/3.0 * gamma * dCoef *  (
                            bc(iStart * grid.dx, jStart * grid.dy, 0.0, simData.currTime - simData.dt) +
                            bc(iStart * grid.dx, jStart * grid.dy, 0.0, simData.currTime) );

                    rhsC.back() = bc(iStart * grid.dx, jStart * grid.dy, (kStart + 0.5 + matA.getSize() - 1) * grid.dz, simData.currTime);
                    break;
                default:
                    break;
            }   
            
           

            break;
        }

        case BoundaryType::Tangent: {
            diag.front() = 1.0;
            double k = simData.k(fieldComponent).valueWithOffset(
                    iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1);
            double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
            double gamma = simData.dt * simData.nu * 0.5 / beta;
            diag.back() = 1 + 4.0 * gamma * dCoef;
            subdiag.back() = -4.0 / 3.0  * gamma * dCoef;


            std::function<double(double t, double x, double y, double z)> bc;
            double xOff = 0.0, yOff = 0.0, zOff = 0.0;
            switch (fieldComponent) {
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
            rhsC.front() = bc((iStart + xOff) * grid.dx, (jStart + yOff) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime);
            

            switch (derivativeDirection) {
                case Axis::X:
                    rhsC.back() =
                            xi(fieldComponent)
                                    .valueWithOffset(iStart, jStart, kStart,
                                                     derivativeDirection, matA.getSize() - 1) 

                            + 4.0 /3.0 * gamma * dCoef *
                            eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 2) 

                            - 4.0 * gamma * dCoef *
                            eta(fieldComponent)
                                    .valueWithOffset(iStart, jStart, kStart,
                                                     derivativeDirection, matA.getSize() - 1) 
                            + 8.0/3.0 * gamma * dCoef *  (
                            bc((iStart + matA.getSize() -1 + 0.5) * grid.dx,
                               (jStart + yOff) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime - simData.dt) +

                            bc((iStart + matA.getSize() -1 + 0.5) * grid.dx,
                               (jStart + yOff) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime) );
                    break;

                case Axis::Y:
                    rhsC.back() =
                            xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                               matA.getSize() - 1) +
                            + 4.0 / 3.0 * gamma * dCoef *
                            eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 2) +
                            - 4.0 * gamma * dCoef *
                            eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 1) +


                            + 8.0 / 3.0 * gamma * dCoef * 
                            (bc((iStart + xOff) * grid.dx,
                                                      (jStart + matA.getSize() -1 + 0.5) * grid.dy,
                                                      (kStart + zOff) * grid.dz, simData.currTime - simData.dt) +

                            bc((iStart + xOff) * grid.dx,
                               (jStart + matA.getSize() -1 + 0.5) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime) );
                    break;

                case Axis::Z:
                    rhsC.back() =
                            xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                               matA.getSize() - 1) 
                            +4.0 / 3.0 * gamma * dCoef *
                            eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 2) 
                            -4.0 * gamma * dCoef *
                            eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 1) 

                            + 8.0/3.0 * gamma * dCoef * (
                            bc((iStart + xOff) * grid.dx,
                                                      (jStart + yOff) * grid.dy,
                                                      (kStart + matA.getSize() -1 + 0.5) * grid.dz, simData.currTime - simData.dt) +

                            bc((iStart + xOff) * grid.dx, (jStart + yOff) * grid.dy,
                               (kStart + matA.getSize() -1 + 0.5) * grid.dz, simData.currTime) );
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
    // Solves Ax = f using Thomas algorithm.
    const int size = matA.getSize();
    
    // Get diagonal vectors: 'a' by reference (unmodified), 'b' and 'c' as working copies.
    const std::vector<double>& a_ref = matA.getDiag(-1); // Reference (no copy needed)
    std::vector<double> b_work = matA.getDiag(0);        // Working copy for denominators
    std::vector<double> c_work = matA.getDiag(1);        // Working copy for c' coefficients

    // Forward Sweep (Elimination)
    if (std::abs(b_work[0]) < 1e-18) 
        { throw std::runtime_error("LinearSys::ThomaSolver failed: Zero pivot detected at the start of elimination."); }
    
    c_work[0] /= b_work[0];
    rhsC[0] /= b_work[0];

    for (unsigned int i = 1; i < size; ++i) {
        double denominator = b_work[i] - a_ref[i - 1] * c_work[i - 1];
        if (std::abs(denominator) < 1e-18) { /* Handle zero pivot */ }
        
        rhsC[i] = (rhsC[i] - a_ref[i - 1] * rhsC[i - 1]) / denominator;
        
        if (i < size - 1) { 
            c_work[i] = c_work[i] / denominator;
        }
    }

    // Backward Sweep (Substitution)
    unknownX[size - 1] = rhsC[size - 1];
    
    for (int i = size - 2; i >= 0; --i) {
        unknownX[i] = rhsC[i] - c_work[i] * unknownX[i + 1];
    }
}
