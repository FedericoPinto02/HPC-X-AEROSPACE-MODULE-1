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
    PressureStep(SimulationData& simData);
    
    

    /**
     * @brief Run pressure step.
     */
    void run();

private:

     Field psi, phi, pcr, divU;   
      // pcr stands for pressure corrector, it's the second phi
      
    /**
     * @brief Construct temporary fields to proceed in computations
     * @param Grid pointer
     */
    void initializeWorkspaceFields(std::shared_ptr<const Grid> gridPtr);



    SimulationData& data_;



};