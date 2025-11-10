#ifndef NSBSOLVER_DERIVATIVES_HPP
#define NSBSOLVER_DERIVATIVES_HPP

#include <vector>

#include "core/Fields.hpp"

/**
 * @brief Handler for computing spatial derivatives of scalar fields.
 */
class Derivatives {
public:
    void computeGradient(const Field &field, VectorField &gradient) const;

    void computeDx(const Field &field, Field &dx) const;

    void computeDy(const Field &field, Field &dy) const;

    void computeDz(const Field &field, Field &dz) const;

    void computeDxDiv(const Field &field, Field &dx) const;

    void computeDyDiv(const Field &field, Field &dy) const;

    void computeDzDiv(const Field &field, Field &dz) const;

    void computeDivergence(const VectorField &field, Field &divergence) const;

    void computeHessianDiag(const Field &field, VectorField &hessianDiag) const;

    void computeDxx(const Field &field, Field &dxx) const;

    void computeDyy(const Field &field, Field &dyy) const;

    void computeDzz(const Field &field, Field &dzz) const;
};

#endif // NSBSOLVER_DERIVATIVES_HPP