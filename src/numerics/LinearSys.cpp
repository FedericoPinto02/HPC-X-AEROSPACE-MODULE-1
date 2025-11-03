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
    // Assign RHS homogeneous b.c. 
    rhsC.front() = 0.0;
    rhsC.back() = 0.0;


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
    subdiag.back() *= 2;
    std::fill(supdiag.begin(), supdiag.end(), -d);
    supdiag.front() *= 2;
    std::fill(diag.begin(), diag.end(), 1 + 2 * d);
}

void LinearSys::fillSystemVelocity(
    const Field &porosity, const VectorField &eta,
    const VectorField &xi, const VectorField &uBoundNew, const VectorField &uBoundOld,
    const Axis fieldComponent, const Axis derivativeDirection,
    const size_t iStart, const size_t jStart, const size_t kStart, 
    const double nu, const double dt) {

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
            "Dimension mismatch: rhsC must be size n.");
    }
    
    // Need to complete RHS, once onto the switch case!

    std::shared_ptr<const Grid> grid = eta.getGrid();
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
        double k = porosity.valueWithOffset(iStart, jStart, kStart,
                                                derivativeDirection, i);
        double beta = 1 + (dt * nu * 0.5 / k);
        double gamma = dt * nu * 0.5 / beta;
        subdiag[i - 1] = -gamma * dCoef;
        supdiag[i] = -gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;
    }

    switch (boundaryType) {
    case BoundaryType::Normal:

    {
        switch (fieldComponent) {
        case Axis::X: // it's a dirichlet on u component
            rhsC.front() =
                uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart,
                                                   Axis::X, 0) -
                grid->dx / 2 *
                    ((uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Y, 1) -
                      uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Y, -1)) /
                         grid->dy +
                     (uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Z, 1) -
                      uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Z, -1)) /
                         grid->dz);
            break;
        case Axis::Y:
            rhsC.front() =
                uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart,
                                                   Axis::Y, 0) -
                grid->dy / 2 *
                    ((uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::X, 1) -
                      uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::X, -1)) /
                         grid->dx +
                     (uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Z, 1) -
                      uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Z, -1)) /
                         grid->dz);
            break;
        case Axis::Z:
            rhsC.front() =
                uBoundNew(Axis::Z).valueWithOffset(iStart, jStart, kStart,
                                                   Axis::Z, 0) -
                grid->dz / 2 *
                    ((uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Y, 1) -
                      uBoundNew(Axis::Y).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::Y, -1)) /
                         grid->dy +
                     (uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::X, 1) -
                      uBoundNew(Axis::X).valueWithOffset(iStart, jStart, kStart,
                                                         Axis::X, -1)) /
                         grid->dx);
            break;
        default:
            break;
        }
        
        diag.back() = 1.0;  // the other is aready initialized to zero
        diag.front() = 1.0;

        rhsC.back() = uBoundNew(fieldComponent)
                          .valueWithOffset(iStart, jStart, kStart,
                                           derivativeDirection, matA.getSize());
        break;

    case BoundaryType::Tangent: {
        diag.front() = 1.0;
        double k = porosity.valueWithOffset(
            iStart, jStart, kStart, derivativeDirection, matA.getSize());
        double beta = 1 + (dt * nu * 0.5 / k);
        double gamma = dt * nu * 0.5 / beta;
        diag.back() = 1 + 3 * gamma * dCoef;
        subdiag.back() = -gamma * dCoef;

        rhsC.front() = uBoundNew(fieldComponent)
                           .valueWithOffset(iStart, jStart, kStart,
                                            derivativeDirection, 0);

        rhsC.back() =
            2 * gamma * dCoef *
                uBoundNew(fieldComponent)
                    .valueWithOffset(iStart, jStart, kStart,
                                     derivativeDirection, matA.getSize()) +
            2 * gamma * dCoef *
                uBoundOld(fieldComponent)
                    .valueWithOffset(iStart, jStart, kStart,
                                     derivativeDirection, matA.getSize()) +
            gamma * dCoef *
                eta(fieldComponent)
                    .valueWithOffset(iStart, jStart, kStart,
                                     derivativeDirection, matA.getSize() - 1) +
            -3 * gamma * dCoef *
                eta(fieldComponent)
                    .valueWithOffset(iStart, jStart, kStart,
                                     derivativeDirection, matA.getSize()) +
            gamma * dCoef *
                xi(fieldComponent)
                    .valueWithOffset(iStart, jStart, kStart,
                                     derivativeDirection, matA.getSize());
        break;
    }
    default:
        break;
    }
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
