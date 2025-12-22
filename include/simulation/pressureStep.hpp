#pragma once

#include <vector>
#include <memory>

#include "core/HaloHandler.hpp"
#include "core/TridiagMat.hpp"
#include "numerics/derivatives.hpp"
#include "numerics/SchurSolver.hpp"
#include "simulation/initializer.hpp"
#include "simulation/SimulationContext.hpp"


/**
 * @brief Handles all pressure step manipulation.
 * This class does not own data but regulates the workflow.
 */
class PressureStep {
    friend class PressureStepTest;

    friend class PressureStepRobustnessTest;

public:

    /**
     * @brief Constructor that initializes the <code>PressureStep</code> with MPI environment and simulation data.
     * @param mpi the MPI environment
     * @param simData the whole simulation data
     */
    PressureStep(MpiEnv &mpi, SimulationData &simData);

    /// Setup fields, linear systems and solvers.
    void setup();

    /// Run the pressure step: pressure correction computation by splitting direction and pressure update.
    void run();

private:
    // --- Environment and helpers -------------------------------------------------------------------------------------
    const MpiEnv &mpi;
    HaloHandler haloHandler;
    Derivatives deriv;

    // --- Phyisics data (pcr stands for pressure corrector, it's the second phi) --------------------------------------
    SimulationData &data_;
    Field psi, phi, pcr, divU;

    // --- Linear system: O(N) memory overhead, O(1) time setup complexity ---------------------------------------------
    std::unique_ptr<SchurSolver> solver_x, solver_y, solver_z;  // solvers for each direction
    std::vector<double> rhs, solution;                          // scratch vectors for linear system solving

    /**
     * @brief Assemble the local linear system matrix for pressure-like variables along a given direction.
     * @param grid the pointer to the grid information
     * @param direction the sweep direction i.e., the direction of the line to be solved
     * @param matrix the tridiagonal matrix to be assembled
     */
    void assembleLocalMatrix(const GridPtr &grid, const Axis direction, TridiagMat &matrix);
};