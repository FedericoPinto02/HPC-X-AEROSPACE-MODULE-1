#ifndef NSBSOLVER_DERIVATIVES_HPP
#define NSBSOLVER_DERIVATIVES_HPP

#include <vector>

#include "core/Fields.hpp"

/**
 * @brief Handler for computing spatial derivatives of scalar fields.
 */
class Derivatives {
public:
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
    void computeDx(const Field &field, Field &dx) const;

    /**
     * @brief Compute the derivative of a scalar field in the y-direction (forward difference).
     * @param field the input scalar field
     * @param dy the output field to store the derivative in y-direction
     */
    void computeDy(const Field &field, Field &dy) const;

    /**
     * @brief Compute the derivative of a scalar field in the z-direction (forward difference).
     * @param field the input scalar field
     * @param dz the output field to store the derivative in z-direction
     */
    void computeDz(const Field &field, Field &dz) const;

    /**
     * @brief Compute the derivative of a scalar field in the x-direction (backward difference).
     * @param field the input vector field
     * @param dx the output field to store the derivative in x-direction
     */
    void computeDxDiv(const Field &field, Field &dx) const;

    /**
     * @brief Compute the derivative of a scalar field in the y-direction (backward difference).
     * @param field the input vector field
     * @param dy the output field to store the derivative in y-direction
     */
    void computeDyDiv(const Field &field, Field &dy) const;

    /**
     * @brief Compute the derivative of a scalar field in the z-direction (backward difference).
     * @param field the input vector field
     * @param dz the output field to store the derivative in z-direction
     */
    void computeDzDiv(const Field &field, Field &dz) const;

    /**
     * @brief Compute the divergence of a vector field (backward difference).
     * @param field the input vector field
     * @param divergence the output scalar field to store the divergence
     */
    void computeDivergence(const VectorField &field, Field &divergence) const;

    /**
     * @brief Compute the diagonal of the Hessian matrix of a scalar field.
     * @param field the input scalar field
     * @param hessianDiag the output vector field to store the Hessian diagonal components
     */
    void computeHessianDiag(const Field &field, VectorField &hessianDiag) const;

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

    /**
     * @brief Compute the local second derivative ∂²f/∂x² at a single grid point.
     *
     * Evaluates the second–order central finite difference at coordinates (i,j,k):
     * ( f(i+1,j,k) + f(i-1,j,k) − 2 f(i,j,k) ) / dx².
     * Performs minimal boundary checking: if @p i is on the domain boundary, returns 0.0.
     *
     * @param f The input scalar field.
     * @param i Grid index in the x–direction.
     * @param j Grid index in the y–direction.
     * @param k Grid index in the z–direction.
     * @return The second derivative ∂²f/∂x² at (i,j,k).
     */
    double Dxx_local(const Field& f, size_t i, size_t j, size_t k) const;

    /**
     * @brief Compute the local second derivative ∂²f/∂y² at a single grid point.
     *
     * Evaluates the second–order central finite difference at coordinates (i,j,k):
     * ( f(i,j+1,k) + f(i,j-1,k) − 2 f(i,j,k) ) / dy².
     * Performs minimal boundary checking: if @p j is on the domain boundary, returns 0.0.
     *
     * @param f The input scalar field.
     * @param i Grid index in the x–direction.
     * @param j Grid index in the y–direction.
     * @param k Grid index in the z–direction.
     * @return The second derivative ∂²f/∂y² at (i,j,k).
     */
    double Dyy_local(const Field& f, size_t i, size_t j, size_t k) const;

    /**
     * @brief Compute the local second derivative ∂²f/∂z² at a single grid point.
     *
     * Evaluates the second–order central finite difference at coordinates (i,j,k):
     * ( f(i,j,k+1) + f(i,j,k-1) − 2 f(i,j,k) ) / dz².
     * Performs minimal boundary checking: if @p k is on the domain boundary, returns 0.0.
     *
     * @param f The input scalar field.
     * @param i Grid index in the x–direction.
     * @param j Grid index in the y–direction.
     * @param k Grid index in the z–direction.
     * @return The second derivative ∂²f/∂z² at (i,j,k).
     */
    double Dzz_local(const Field& f, size_t i, size_t j, size_t k) const;

};

#endif // NSBSOLVER_DERIVATIVES_HPP