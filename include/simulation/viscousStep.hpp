# pragma once

#include <vector>
#include <memory>
#include <simulation/initializer.hpp>
#include <numerics/LinearSys.hpp>
#include <numerics/derivatives.hpp>
#include <simulation/SimulationContext.hpp>


/**
 * @brief Handles all viscous step maniplation.
 * This class does not own data but regulates the workflow.
 */
class ViscousStep
{
public:

    /**
     * @brief Constructor
     * @param contex Contains all simulation data
     */
    ViscousStep(SimulationContext& ctx);
    
    

    /**
     * @brief Run viscous step.
     */
    void run();

private:

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

    SimulationContext& context_;



};