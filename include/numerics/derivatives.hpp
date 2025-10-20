#ifndef NSBSOLVER_DERIVATIVES_HPP
#define NSBSOLVER_DERIVATIVES_HPP

#include <vector>

#include "core/Fields.hpp"

/**
 * @brief Handler for computing spatial derivatives of scalar fields.
 */
class Derivatives {
public:
    void computeGradient(const Field &field, VectorField &gradient);

    void computeDx(const Field &field, Field &dx);

    void computeDy(const Field &field, Field &dy);

    void computeDz(const Field &field, Field &dz);

    void computeDivergence(const Field &field, Field &divergence);

    void computeHessianDiag(const Field &field, VectorField &hessianDiag);

    void computeDxx(const Field &field, Field &dxx);

    void computeDyy(const Field &field, Field &dyy);

    void computeDzz(const Field &field, Field &dzz);
};

#endif // NSBSOLVER_DERIVATIVES_HPP