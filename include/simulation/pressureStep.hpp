# pragma once

#include <vector>
#include <memory>
#include <simulation/initializer.hpp>
#include <numerics/LinearSys.hpp>
#include <numerics/derivatives.hpp>
#include <simulation/SimulationContext.hpp>


/**
 * @brief Handles all pressure step manipulation.
 * This class does not own data but regulates the workflow.
 */
class PressureStep
{
    friend class PressureStepTest;
    friend class PressureStepRobustnessTest;
public:

    /**
     * @brief Constructor.
     * @param simData the simulation data
     * @param parallel the parallelization settings
     */
    PressureStep(SimulationData& simData, ParallelizationSettings& parallel);

    /// Run the pressure step: pressure correction computation by splitting direction and pressure update.
    void run();

private:
    SimulationData& data_;
    ParallelizationSettings parallel_;

    Field psi, phi, pcr, divU;
    // pcr stands for pressure corrector, it's the second phi

    /**
     * @brief Wrapper that chooses whether to use Thomas (P=1) or Schur (P>1).
     */
    std::vector<double> solveSystem(LinearSys& sys);

};