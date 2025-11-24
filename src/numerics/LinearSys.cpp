#include <algorithm>
#include <cmath>
#include <core/Fields.hpp>
#include <numerics/LinearSys.hpp>
#include <numerics/derivatives.hpp>
#include <vector>

LinearSys::LinearSys(int n, BoundaryType boundaryType)
        : boundaryType(boundaryType), matA(n) {
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

void LinearSys::fillSystemPressure(const Field &phi, const Axis direction) {

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
                "Dimension mismatch: rhsIncomplete must be size n.");
    }


    const Grid &grid = phi.getGrid();
    double d = 0;
    switch (direction) {
        case Axis::X:
            d = 1 / grid.dx / grid.dx;
            break;
        case Axis::Y:
            d = 1 / grid.dy / grid.dy;
            break;
        case Axis::Z:
            d = 1 / grid.dz / grid.dz;
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
        const SimulationData &simData, const VectorField &xi,
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

    const Grid &grid = simData.eta.getGrid();
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
        double k = simData.k.valueWithOffset(iStart, jStart, kStart,
                                             derivativeDirection, i);
        double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
        double gamma = simData.dt * simData.nu * 0.5 / beta;
        subdiag[i - 1] = -gamma * dCoef;
        supdiag[i] = -gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;
    }

    switch (boundaryType) {
        case BoundaryType::Normal: {
            switch (fieldComponent) {
                case Axis::X:
                    rhsC.front() =
                            simData.bcu(0.0, jStart * grid.dy, kStart * grid.dz, simData.currTime) -
                            grid.dx / 2 *
                            ((simData.bcv(0.0, (jStart + 0.5) * grid.dy, kStart * grid.dz, simData.currTime) -
                              simData.bcv(0.0, (jStart - 0.5) * grid.dy, kStart * grid.dz, simData.currTime)) /
                             grid.dy +
                             (simData.bcw(0.0, jStart * grid.dy, (kStart + 0.5) * grid.dz, simData.currTime) -
                              simData.bcw(0.0, jStart * grid.dy, (kStart - 0.5) * grid.dz, simData.currTime)) /
                             grid.dz);
                    break;

                case Axis::Y:
                    rhsC.front() =
                            simData.bcv(iStart * grid.dx, 0.0, kStart * grid.dz, simData.currTime) -
                            grid.dy / 2 *
                            ((simData.bcu((iStart + 0.5) * grid.dx, 0.0, kStart * grid.dz, simData.currTime) -
                              simData.bcu((iStart - 0.5) * grid.dx, 0.0, kStart * grid.dz, simData.currTime)) /
                             grid.dx +
                             (simData.bcw(iStart * grid.dx, 0.0, (kStart + 0.5) * grid.dz, simData.currTime) -
                              simData.bcw(iStart * grid.dx, 0.0, (kStart - 0.5) * grid.dz, simData.currTime)) /
                             grid.dz);
                    break;

                case Axis::Z:
                    rhsC.front() =
                            simData.bcw(iStart * grid.dx, jStart * grid.dy, 0.0, simData.currTime) -
                            grid.dz / 2 *
                            ((simData.bcu((iStart + 0.5) * grid.dx, jStart * grid.dy, 0.0, simData.currTime) -
                              simData.bcu((iStart - 0.5) * grid.dx, jStart * grid.dy, 0.0, simData.currTime)) /
                             grid.dx +
                             (simData.bcv(iStart * grid.dx, (jStart + 0.5) * grid.dy, 0.0, simData.currTime) -
                              simData.bcv(iStart * grid.dx, (jStart - 0.5) * grid.dy, 0.0, simData.currTime)) /
                             grid.dy);
                    break;

                default:
                    break;
            }

            diag.back() = 1.0; // the other is aready initialized to zero
            diag.front() = 1.0;
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
            
            switch (derivativeDirection) {
                case Axis::X:
                    rhsC.back() = bc((iStart + xOff + matA.getSize() - 1) * grid.dx, (jStart + yOff) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime);
                    break;

                case Axis::Y:
                     rhsC.back() = bc((iStart + xOff) * grid.dx, (jStart + yOff + matA.getSize() - 1) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime);
                    break;

                case Axis::Z:
                     rhsC.back() = bc((iStart + xOff) * grid.dx, (jStart + yOff) * grid.dy, (kStart + zOff + matA.getSize() - 1) * grid.dz, simData.currTime);
                    break;
                default:
                    break;
            }   

            break;
        }

        case BoundaryType::Tangent: {
            diag.front() = 1.0;
            double k = simData.k.valueWithOffset(
                    iStart, jStart, kStart, derivativeDirection, matA.getSize() - 1);
            double beta = 1 + (simData.dt * simData.nu * 0.5 / k);
            double gamma = simData.dt * simData.nu * 0.5 / beta;
            diag.back() = 1 + 3 * gamma * dCoef;
            subdiag.back() = -gamma * dCoef;


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
                                                     derivativeDirection, matA.getSize() - 1) +

                            -gamma * dCoef *
                            simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 2) +

                            3.0 * gamma * dCoef *
                            simData.eta(fieldComponent)
                                    .valueWithOffset(iStart, jStart, kStart,
                                                     derivativeDirection, matA.getSize() - 1) +
                            -2.0 * gamma * dCoef *  // uOld
                            bc((iStart + matA.getSize() - 1) * grid.dx,
                               (jStart + yOff) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime - simData.dt) +

                            2.0 * gamma * dCoef *
                            bc((iStart + matA.getSize() - 1) * grid.dx,
                               (jStart + yOff) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime);
                    break;

                case Axis::Y:
                    rhsC.back() =
                            xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                               matA.getSize() - 1) +
                            -gamma * dCoef *
                            simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 2) +
                            3.0 * gamma * dCoef *
                            simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 1) +


                            -2.0 * gamma * dCoef * bc((iStart + xOff) * grid.dx,
                                                      (jStart + matA.getSize() - 1) * grid.dy,
                                                      (kStart + zOff) * grid.dz, simData.currTime - simData.dt) +

                            2.0 * gamma * dCoef *
                            bc((iStart + xOff) * grid.dx,
                               (jStart + matA.getSize() - 1) * grid.dy, (kStart + zOff) * grid.dz, simData.currTime);
                    break;

                case Axis::Z:
                    rhsC.back() =
                            xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                               matA.getSize() - 1) +
                            -gamma * dCoef *
                            simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 2) +
                            3.0 * gamma * dCoef *
                            simData.eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection,
                                                                        matA.getSize() - 1) +

                            -2.0 * gamma * dCoef * bc((iStart + xOff) * grid.dx,
                                                      (jStart + yOff) * grid.dy,
                                                      (kStart + matA.getSize() - 1) * grid.dz, simData.currTime - simData.dt) +

                            2.0 * gamma * dCoef *
                            bc((iStart + xOff) * grid.dx, (jStart + yOff) * grid.dy,
                               (kStart + matA.getSize() - 1) * grid.dz, simData.currTime);
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
