#ifndef NSBSOLVER_DERIVATIVES_HPP
#define NSBSOLVER_DERIVATIVES_HPP

#include <vector>

#include "core/Fields.hpp"

/**
 * @brief Handler for computing spatial derivatives of scalar fields.
 */
class Derivatives {
public:
    void computeGradient(const Field &field, VectorField &gradient) {
        computeDx(field, gradient.x());
        computeDy(field, gradient.y());
        computeDz(field, gradient.z());
        // todo - may be optimized in terms of locality: tiling ??
    }

    void computeDx(const Field &field, Field &dx);

    void computeDy(const Field &field, Field &dy);

    void computeDz(const Field &field, Field &dz);

    void computeDivergence(const Field &field, Field &divergence) {
        Field tmp;
        tmp.setup(field.getGrid(),
                  std::vector<Field::Scalar>(field.getGrid()->size(), 0.0));

        divergence.reset();
        computeDx(field, tmp);
        divergence.add(tmp);
        computeDy(field, tmp);
        divergence.add(tmp);
        computeDz(field, tmp);
        divergence.add(tmp);
    }

    void computeHessianDiag(const Field &field, VectorField &hessianDiag) {
        computeDx(field, hessianDiag.x());
        computeDy(field, hessianDiag.y());
        computeDz(field, hessianDiag.z());
        // todo - may be optimized in terms of locality: tiling ??
    }

    void computeDxx(const Field &field, Field &dxx);

    void computeDyy(const Field &field, Field &dyy);

    void computeDzz(const Field &field, Field &dzz);
};

#endif // NSBSOLVER_DERIVATIVES_HPP