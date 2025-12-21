# pragma once

#include <memory>
#include <vector>

#include "core/HaloHandler.hpp"
#include "core/TridiagMat.hpp"
#include "numerics/derivatives.hpp"
#include "simulation/initializer.hpp"
#include "numerics/SchurSolver.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Handles all viscous step maniplation.
 * This class does not own data but regulates the workflow.
 */
class ViscousStep
{
    friend class ViscousStepTest;

    friend class ViscousStepRobustnessTest;
public:

    /**
     * @brief Constructor that initializes the <code>ViscousStep</code> with MPI environment and simulation data.
     * @param mpi the MPI environment
     * @param simData the whole simulation data
     */
    ViscousStep(MpiEnv &mpi, SimulationData &simData);

    /// Run viscous step.
    void run();

private:
    MpiEnv &mpi;
    SimulationData& data_;

    VectorField g, gradP, dxxEta, dyyZeta, dzzU, xi;

    /// Compute the g term (necessary for the xi term).
    void computeG();

    /// Compute the xi term (necessary for the x-sweep).
    void computeXi();

    /// Solves the viscous step by ADI method (x-, y- and z-sweeps) and updates the velocity field.
    void closeViscousStep();

    /**
     * @brief Assemble the local linear system for velocity-like variables, for a given field component
     * and derivative direction.
     * @param simData the whole simulation data
     * @param eta the unknown vector field at the previous time step
     * @param xi the known vector field from a previous sweep at the current time step
     * @param fieldComponent the component of the vector field to solve for
     * @param derivativeDirection the sweep direction i.e., the direction of the second derivative
     * @param iStart the starting index in the i-direction
     * @param jStart the starting index in the j-direction
     * @param kStart the starting index in the k-direction
     * @param matrix the tridiagonal matrix to be assembled
     * @param rhs the right-hand side vector to assemble boundary conditions into
     */
    void assembleLocalSystem(
            const SimulationData &simData, const VectorField &eta, const VectorField &xi,
            const Axis fieldComponent, const Axis derivativeDirection,
            const size_t iStart, const size_t jStart, const size_t kStart,
            TridiagMat &matrix, std::vector<double> &rhs
    );
};