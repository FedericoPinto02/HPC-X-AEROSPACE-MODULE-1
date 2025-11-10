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


void NSBSolver::solve() {

    // ----- 1. Read input -----
    InputReader reader;
    InputData input = reader.read("../data/config.json");

    std::cout << "[OK] Input file read successfully.\n";
    LogWriter logger(input.logging);
    logger.printInputSummary(input);

    // ----- 2. Simulation Setup -----
    std::cout << "Initializing simulation setup...\n";
    Initializer init;
    SimulationData simData = init.setup(input);

    std::cout << "[OK] SimulationData created successfully.\n\n";
    logger.printRuntimeSummary(simData);

    // Initialize steps & VTK writer
    OutputSettings output = input.output;
    LoggingSettings logging = input.logging;
    ViscousStep viscousStep{simData};
    PressureStep pressureStep{simData};
    VTKWriter vtkWriter(simData.Nx, simData.Ny, simData.Nz,
              simData.dx, simData.dy, simData.dz);


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
        viscousStep.run();
        pressureStep.run();
        auto end = std::chrono::high_resolution_clock::now();    // stop timer

        std::chrono::duration<double> elapsed = end - start;
        if (simData.currStep % logging.loggingFrequency == 0) {
            logger.printLoopProgress(simData, elapsed.count());
        }
        
        // --- ON-THE-FLY OUTPUT LOGIC ---

        // 1. Check if we should write output this step
        if (output.enabled && (simData.currStep % output.outputFrequency == 0))
        {
            // 2. Prepare arguments
            std::string baseprefix = output.dir + "/" + output.baseFilename;
            int step = static_cast<int>(simData.currStep);
            std::string title = "Simulation Data";

            // 3. Create shared_ptrs "on the fly" (makes copies)
            auto pressure_ptr = std::make_shared<Field>(simData.p);
            auto velocity_ptr = std::make_shared<VectorField>(simData.u);

            // 4. Call the 5-argument writer function
            vtkWriter.write_timestep(baseprefix, 
                                     step, 
                                     pressure_ptr, 
                                     velocity_ptr, 
                                     title);
            
            std::cout << "VTK file written for step " << step << "\n";
        }
        // --- END OUTPUT LOGIC ---

    }
}