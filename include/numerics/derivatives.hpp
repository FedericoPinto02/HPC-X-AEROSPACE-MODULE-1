#pragma once

#include <vector>

#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @class Derivatives
 * @brief Computational engine for spatial derivatives on a staggered or collocated grid.
 * * This class provides methods to compute first and second-order spatial derivatives 
 * using finite difference stencils. It handles:
 * - **Forward Differences**
 * - **Backward Differences**
 * - **Centered Differences**
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
     * @param bcu Boundary condition function for the U component.
     * @param time Current simulation time (for time-dependent BCs).
     */
    void computeDx_bwd(const Field &field, Field &dx, Func &bcu, double time) const;

    /**
     * @brief Compute the derivative of a scalar field in the y-direction (backward difference).
     * @param field the input vector field
     * @param dy the output field to store the derivative in y-direction
     * @param bcv Boundary condition function for the V component.
     * @param time Current simulation time.
     */
    void computeDy_bwd(const Field &field, Field &dy, Func &bcv, double time) const;

    /**
     * @brief Compute the derivative of a scalar field in the z-direction (backward difference).
     * @param field the input vector field
     * @param dz the output field to store the derivative in z-direction
     * @param bcw Boundary condition function for the W component.
     * @param time Current simulation time.
     */
    void computeDz_bwd(const Field &field, Field &dz, Func &bcw, double time) const;

    /**
     * @brief Compute the divergence of a vector field (backward difference).
     * @param field the input vector field
     * @param divergence the output scalar field to store the divergence
     * @param bcu, bcv, bcw Boundary functions for each velocity component.
     * @param time Current simulation time.
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
     * @param bcu, bcv, bcw Boundary condition functions.
     * @param time Current time. 
     */
    void computeDxx(const VectorField &field, VectorField &dxx, const Func &bcu, const Func &bcv, const Func &bcw, const double &time) const;

    /**
     * @brief Compute the second derivative of a vector field in the y-direction.
     * @param field the input vector field
     * @param dyy the output vector field to store the second derivative in y-direction
     * @param bcu, bcv, bcw Boundary condition functions.
     * @param time Current time. 
     */
    void computeDyy(const VectorField &field, VectorField &dyy, const Func &bcu, const Func &bcv, const Func &bcw, const double &time) const;

    /**
     * @brief Compute the second derivative of a vector field in the z-direction.
     * @param field the input vector field
     * @param dzz the output vector field to store the second derivative in z-direction
     * @param bcu, bcv, bcw Boundary condition functions.
     * @param time Current time. 
     */
    void computeDzz(const VectorField &field, VectorField &dzz, const Func &bcu, const Func &bcv, const Func &bcw, const double &time) const;

    /**
     * @brief Compute the second derivative of a scalar field in the x-direction.
     * @param field the input scalar field
     * @param dxx the output field to store the second derivative in x-direction
     * @param bcu, bcv, bcw Boundary condition functions.
     * @param time Current time.
     * @param axis Specifies which component of the velocity/field the BC refers to.
     */
    void computeDxx(const Field &field, Field &dxx, const Func &bc, const double &time, const Axis axis) const;

    /**
     * @brief Compute the second derivative of a scalar field in the y-direction.
     * @param field the input scalar field
     * @param dyy the output field to store the second derivative in y-direction
     * @param bcu, bcv, bcw Boundary condition functions.
     * @param time Current time.
     * @param axis Specifies which component of the velocity/field the BC refers to.
     */
    void computeDyy(const Field &field, Field &dyy, const Func &bc, const double &time, const Axis axis) const;

    /**
     * @brief Compute the second derivative of a scalar field in the z-direction.
     * @param field the input scalar field
     * @param dzz the output field to store the second derivative in z-direction
     * @param bcu, bcv, bcw Boundary condition functions.
     * @param time Current time.
     * @param axis Specifies which component of the velocity/field the BC refers to.
     */
    void computeDzz(const Field &field, Field &dzz, const Func &bc, const double &time, const Axis axis) const;
};
