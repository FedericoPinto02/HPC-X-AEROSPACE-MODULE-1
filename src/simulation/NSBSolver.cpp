#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

// --- Include headers ---
#include "simulation/NSBSolver.hpp"

#include "core/Fields.hpp"
#include "io/inputReader.hpp"
#include "io/VTKWriter.hpp"
#include "io/logWriter.hpp"
#include "numerics/derivatives.hpp"
#include "simulation/pressureStep.hpp"
#include "simulation/SimulationContext.hpp"
#include "simulation/viscousStep.hpp"
#include "simulation/initializer.hpp"

using namespace MuParserXAdapter; 

NSBSolver::NSBSolver(const std::string& configFile) 
    : configFile(configFile) {}

    
void NSBSolver::setup() {

    // Read input
    InputReader reader;
    input = reader.read(configFile);

    // Logger
    logger = std::make_unique<LogWriter>(input.logging);
    
    // Init SimulationData
    Initializer init;
    simData = init.setup(input);

    ParallelizationSettings parallel;
    parallel.schurDomains = input.parallelization.schurDomains;

    // Steps
    viscousStep  = std::make_unique<ViscousStep>(simData, parallel);
    pressureStep = std::make_unique<PressureStep>(simData, parallel);

    // Writer
    vtkWriter = std::make_unique<VTKWriter>(input.output, simData);

    // Store settings
    outputSettings  = input.output;
    loggingSettings = input.logging;

    // Boundary conditions expressions
    simData.bcu = createTimeFunction(input.boundary_conditions.u_expr);
    simData.bcv = createTimeFunction(input.boundary_conditions.v_expr);
    simData.bcw = createTimeFunction(input.boundary_conditions.w_expr);

    // Store initializer
    initializer = std::make_unique<Initializer>(init);

    // Logging
    logger->printSimulationHeader(input, simData);
}


void NSBSolver::solve() {

    if (!initializer) {
        std::cerr << "[ERROR] NSBSolver::setup() was not called before solve().\n";
        std::exit(EXIT_FAILURE);
    }

    auto& init = *initializer;

    logger->printStepHeader();

    // Time integration
    for (unsigned int i = 1; i < simData.totalSteps + 1; i++)
    {
        // 1. Update Time & Fields
        simData.uBoundOld = simData.uBoundNew;
        simData.currStep++;
        simData.currTime += simData.dt; 

        simData.uBoundNew = init.initializeTemporalVectorFieldFromExpr(
            simData.currTime,
            simData.gridPtr,
            input.boundary_conditions.u_expr,
            input.boundary_conditions.v_expr,
            input.boundary_conditions.w_expr);
            
        simData.f = init.initializeTemporalVectorFieldFromExpr(
            simData.currTime,
            simData.gridPtr,
            input.forces.fx_expr,
            input.forces.fy_expr,
            input.forces.fz_expr);

        // 2. Physics Steps (Timed)
        auto start = std::chrono::high_resolution_clock::now(); 
        
        viscousStep->run();
        pressureStep->run();
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        // 3. Output
        bool vtkWritten = vtkWriter->write_timestep_if_needed(
                            simData.currStep, 
                            simData.p, 
                            simData.u);

        // 4. Logging Progress (Tabular)
        logger->printStepProgress(
                simData.currStep,
                simData.currTime,
                simData.dt,
                elapsed.count(),
                vtkWritten
            );
    }
}
