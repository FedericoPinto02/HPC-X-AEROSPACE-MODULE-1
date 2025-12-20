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
     * @brief Constructor
     * @param contex Contains all simulation data
     */
    ViscousStep(MpiEnv &mpi, SimulationData &simData);
    
    /**
     * @brief Run viscous step.
     */
    void run();

private:
    MpiEnv &mpi;
    SimulationData& data_;

    VectorField g, gradP, dxxEta, dyyZeta, dzzU, xi;

    /**
     * @brief Compute G term
     */
    void computeG();

     /**
     * @brief Compute xi term
     */
    void computeXi();

    /**
     * @brief Closes viscous step filling and solving three linear systems.
     */
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