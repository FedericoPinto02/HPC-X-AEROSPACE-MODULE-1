#ifndef NSBSOLVER_DERIVATIVES_HPP
#define NSBSOLVER_DERIVATIVES_HPP

#include <vector>

#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Handler for computing spatial derivatives of scalar fields.
 */
class Derivatives {
public:
    //==================================================================================================================
    //--- First-order forward differences ------------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Compute the gradient of a scalar field (forward difference).
     * @param field the input scalar field
     * @param gradient the output vector field to store the gradient
     */
    void computeGradient(const Field &field, VectorField &gradient) const;

    /**
     * @brief Compute the derivative of a scalar field in the x-direction (forward difference).
     * @param field the input scalar field
     * @param dx the output field to store the derivative in x-direction
     */
    void computeDx_fwd(const Field &field, Field &dx) const;

    /**
     * @brief Compute the derivative of a scalar field in the y-direction (forward difference).
     * @param field the input scalar field
     * @param dy the output field to store the derivative in y-direction
     */
    void computeDy_fwd(const Field &field, Field &dy) const;

    /**
     * @brief Compute the derivative of a scalar field in the z-direction (forward difference).
     * @param field the input scalar field
     * @param dz the output field to store the derivative in z-direction
     */
    void computeDz_fwd(const Field &field, Field &dz) const;

    //==================================================================================================================
    //--- First-order backward differences -----------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Compute the derivative of a scalar field in the x-direction (backward difference).
     * @param field the input vector field
     * @param dx the output field to store the derivative in x-direction
     */
    void computeDx_bwd(const Field &field, Field &dx, Func &bcu, double time) const;

    /**
     * @brief Compute the derivative of a scalar field in the y-direction (backward difference).
     * @param field the input vector field
     * @param dy the output field to store the derivative in y-direction
     */
    void computeDy_bwd(const Field &field, Field &dy, Func &bcv, double time) const;

    /**
     * @brief Compute the derivative of a scalar field in the z-direction (backward difference).
     * @param field the input vector field
     * @param dz the output field to store the derivative in z-direction
     */
    void computeDz_bwd(const Field &field, Field &dz, Func &bcw, double time) const;

    /**
     * @brief Compute the divergence of a vector field (backward difference).
     * @param field the input vector field
     * @param divergence the output scalar field to store the divergence
     */
    void computeDivergence(
            const VectorField &field, Field &divergence,
            Func &bcu, Func &bcv, Func &bcw, double time
    ) const;

    //==================================================================================================================
    //--- Second-order centered differences ----------------------------------------------------------------------------
    //==================================================================================================================
    /**
     * @brief Compute the second derivative of a vector field in the x-direction.
     * @param field the input vector field
     * @param dxx the output vector field to store the second derivative in x-direction
     */
    void computeDxx(const VectorField &field, VectorField &dxx) const;

    /**
     * @brief Compute the second derivative of a vector field in the y-direction.
     * @param field the input vector field
     * @param dyy the output vector field to store the second derivative in y-direction
     */
    void computeDyy(const VectorField &field, VectorField &dyy) const;

    /**
     * @brief Compute the second derivative of a vector field in the z-direction.
     * @param field the input vector field
     * @param dzz the output vector field to store the second derivative in z-direction
     */
    void computeDzz(const VectorField &field, VectorField &dzz) const;

    /**
     * @brief Compute the second derivative of a scalar field in the x-direction.
     * @param field the input scalar field
     * @param dxx the output field to store the second derivative in x-direction
     */
    void computeDxx(const Field &field, Field &dxx) const;

    /**
     * @brief Compute the second derivative of a scalar field in the y-direction.
     * @param field the input scalar field
     * @param dyy the output field to store the second derivative in y-direction
     */
    void computeDyy(const Field &field, Field &dyy) const;

    /**
     * @brief Compute the second derivative of a scalar field in the z-direction.
     * @param field the input scalar field
     * @param dzz the output field to store the second derivative in z-direction
     */
    void computeDzz(const Field &field, Field &dzz) const;
};

#endif // NSBSOLVER_DERIVATIVES_HPP