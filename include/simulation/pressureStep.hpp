# pragma once

#include <vector>
#include <memory>
#include <simulation/initializer.hpp>
#include <numerics/LinearSys.hpp>
#include <numerics/derivatives.hpp>
#include <simulation/SimulationContext.hpp>


/**
 * @brief Handles all pressure step maniplation.
 * This class does not own data but regulates the workflow.
 */
class PressureStep
{
    friend class PressureStepTest;
    friend class PressureStepRobustnessTest;
public:

    /**
     * @brief Constructor
     * @param contex Contains all simulation data
     */
    PressureStep(SimulationData& simData, ParallelizationSettings& parallel);

    /**
     * @brief Run pressure step.
     */
    void run();

private:
    SimulationData& data_;
    ParallelizationSettings parallel_;

    Field psi, phi, pcr, divU;
    // pcr stands for pressure corrector, it's the second phi

    LinearSys linearSys_xLine, linearSys_yLine, linearSys_zLine;


    /**
     * @brief Assembles the tridiagonal matrix for 1D line solver.
     *
     * The tridiagonal stencil is [-1/delta^2, 1 + 2/delta^2, -1/delta^2].
     * @param delta the grid spacing in the line direction
     * @param sys the LinearSys object to assemble the matrix for
     */
    void assembleLineMatrix(double delta, LinearSys& sys);

    /**
     * @brief Wrapper that chooses whether to use Thomas (P=1) or Schur (P>1).
     */
    std::vector<double> solveSystem(LinearSys& sys, BoundaryType bType);

};