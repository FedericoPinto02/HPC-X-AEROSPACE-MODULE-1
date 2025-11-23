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
    friend class ViscousStepTest;
    friend class ViscousStepSolverTest;
    friend class ViscousStepRobustnessTest;
public:

    /**
     * @brief Constructor
     * @param contex Contains all simulation data
     */
    ViscousStep(SimulationData& simData, ParallelizationSettings& parallel);
    
    

    /**
     * @brief Run viscous step.
     */
    void run();

private:

     VectorField g, gradP, dxxEta, dyyZeta, dzzU, xi;    
    /**
     * @brief Construct temporary fields to proceed in computations
     * @param Grid pointer
     */
    void initializeWorkspaceFields(std::shared_ptr<const Grid> gridPtr);

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

    SimulationData& data_;
    ParallelizationSettings parallel_;

    /**
     * @brief Wrapper that chooses whether to use Thomas (P=1) or Schur (P>1).
     */
    std::vector<double> solveSystem(LinearSys& sys, BoundaryType bType);



};