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
using namespace MuParserXAdapter; // / createTimeFunction


NSBSolver::NSBSolver(const std::string& configFile) 
    : configFile(configFile) {}

    
void NSBSolver::setup() {
    // 1. Read input
    InputReader reader;
    input = reader.read(configFile);

    std::cout << "[OK] Input file read successfully.\n";

    // 2. Logger
    logger = std::make_unique<LogWriter>(input.logging);
    logger->printInputSummary(input);

    // 3. Init SimulationData
    std::cout << "Initializing simulation setup...\n";
    Initializer init;
    simData = init.setup(input);

    std::cout << "[OK] SimulationData created successfully.\n\n";
    logger->printRuntimeSummary(simData);

    // 4. Steps
    viscousStep  = std::make_unique<ViscousStep>(simData);
    pressureStep = std::make_unique<PressureStep>(simData);

    // 5. Writer
    vtkWriter = std::make_unique<VTKWriter>(
        simData.Nx, simData.Ny, simData.Nz,
        simData.dx, simData.dy, simData.dz
    );

    // Store output/logging settings
    outputSettings  = input.output;
    loggingSettings = input.logging;

    // ----- 5. Boundary conditions expressions -----
    simData.bcu = createTimeFunction(input.boundary_conditions.u_expr);
    simData.bcv = createTimeFunction(input.boundary_conditions.v_expr);
    simData.bcw = createTimeFunction(input.boundary_conditions.w_expr);

    // Store initializer for temporal fields
    initializer = std::make_unique<Initializer>(init);

}


void NSBSolver::solve() {

    if (!initializer) {
        std::cerr << "[ERROR] NSBSolver::setup() was not called before solve().\n";
        std::exit(EXIT_FAILURE);
    }

    auto& init = *initializer;

    // Time integration
    for (unsigned int i = 1; i < simData.totalSteps + 1; i++)
    {
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
        auto start = std::chrono::high_resolution_clock::now();  // start timer
        viscousStep->run();
        pressureStep->run();
        auto end = std::chrono::high_resolution_clock::now();    // stop timer

        std::chrono::duration<double> elapsed = end - start;
        if (simData.currStep % loggingSettings.loggingFrequency == 0) {
            logger->printLoopProgress(simData, elapsed.count());
        }
        
        // --- ON-THE-FLY OUTPUT LOGIC ---

        // 1. Check if we should write output this step
        if (outputSettings.enabled && (simData.currStep % outputSettings.outputFrequency == 0))
        {
            // 2. Prepare arguments
            std::string baseprefix = outputSettings.dir + "/" + outputSettings.baseFilename;
            int step = static_cast<int>(simData.currStep);
            std::string title = "Simulation Data";

            // 3. Create shared_ptrs "on the fly" (makes copies)
            auto pressure_ptr = std::make_shared<Field>(simData.p);
            auto velocity_ptr = std::make_shared<VectorField>(simData.u);

            // 4. Call the 5-argument writer function
            vtkWriter->write_timestep(baseprefix, 
                                     step, 
                                     pressure_ptr, 
                                     velocity_ptr, 
                                     title);
            
            std::cout << "VTK file written for step " << step << "\n";
        }
        // --- END OUTPUT LOGIC ---

    }
}
