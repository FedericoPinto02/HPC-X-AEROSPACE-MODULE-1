#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>

// --- Include headers ---
#include "simulation/NSBSolver.hpp"

#include "core/Fields.hpp"
#include "io/inputReader.hpp"
#include "io/VTKWriter.hpp"
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

    // ----- 2. Print verbose input summary -----
    // This will only print if "verbose: true" is set in config.json
    if (input.logging.verbose) 
    {
        std::cout << "\n=== INPUT CONFIG SUMMARY ===\n" << std::boolalpha;

        std::cout << "[MESH]\n";
        std::cout << "  - Grid size (Nx, Ny, Nz): " << input.mesh.nx << ", " << input.mesh.ny << ", " << input.mesh.nz << "\n";
        std::cout << "  - Spacing (dx, dy, dz):   " << input.mesh.dx << ", " << input.mesh.dy << ", " << input.mesh.dz << "\n";

        std::cout << "[TIME]\n";
        std::cout << "  - Time step (dt):   " << input.time.dt << "\n";
        std::cout << "  - End time (t_end): " << input.time.t_end << "\n";

        std::cout << "[PHYSICS]\n";
        std::cout << "  - Viscosity (nu):   " << input.physics.nu << "\n";
        std::cout << "  - Permeability (k): \"" << input.physics.k_expr << "\"\n";

        std::cout << "[FORCES]\n";
        std::cout << "  - Force x (fx): \"" << input.forces.fx_expr << "\"\n";
        std::cout << "  - Force y (fy): \"" << input.forces.fy_expr << "\"\n";
        std::cout << "  - Force z (fz): \"" << input.forces.fz_expr << "\"\n";

        std::cout << "[INITIAL CONDITIONS]\n";
        std::cout << "  - Velocity u: \"" << input.initial_conditions.u_expr << "\"\n";
        std::cout << "  - Velocity v: \"" << input.initial_conditions.v_expr << "\"\n";
        std::cout << "  - Velocity w: \"" << input.initial_conditions.w_expr << "\"\n";
        std::cout << "  - Pressure p: \"" << input.initial_conditions.p_expr << "\"\n";

        std::cout << "[BOUNDARY CONDITIONS]\n";
        std::cout << "  - Velocity u: \"" << input.boundary_conditions.u_expr << "\"\n";
        std::cout << "  - Velocity v: \"" << input.boundary_conditions.v_expr << "\"\n";
        std::cout << "  - Velocity w: \"" << input.boundary_conditions.w_expr << "\"\n";

        std::cout << "[OUTPUT]\n";
        std::cout << "  - Enabled:    " << input.output.enabled << "\n";
        std::cout << "  - Directory:  \"" << input.output.dir << "\"\n";
        std::cout << "  - Filename:   \"" << input.output.baseFilename << "\"\n";
        std::cout << "  - Frequency:  " << input.output.frequency << "\n";

        std::cout << "[LOGGING]\n";
        std::cout << "  - Verbose:        " << input.logging.verbose << "\n";
        std::cout << "  - Log to File:    " << input.logging.logToFile << "\n";
        std::cout << "  - Log to Console: " << input.logging.logToConsole << "\n";
        std::cout << "  - Filename:       \"" << input.logging.filename << "\"\n";
        std::cout << "  - Frequency:      " << input.logging.frequency << "\n";
        
        std::cout << "============================\n\n";
    }

    // ----- 3. Simulation Setup -----
    std::cout << "Initializing simulation setup...\n";
    Initializer init;
    SimulationData simData = init.setup(input);
    std::cout << "[OK] SimulationData created successfully.\n\n";

    // ----- 4. Print simulation summary -----
    std::cout << "=== SIMULATION RUNTIME SETTINGS ===\n";
    
    std::cout << "[GRID]\n";
    std::cout << "  - Grid size (Nx, Ny, Nz): " << simData.Nx << " x " << simData.Ny << " x " << simData.Nz << "\n";
    std::cout << "  - Domain (Lx, Ly, Lz):  " << simData.Lx << ", " << simData.Ly << ", " << simData.Lz << "\n";
    std::cout << "  - Spacing (dx, dy, dz): " << simData.dx << ", " << simData.dy << ", " << simData.dz << "\n";

    std::cout << "[TIME]\n";
    std::cout << "  - Time step (dt): " << simData.dt << "\n";
    std::cout << "  - Total time:     " << simData.totalSimTime << "\n";
    std::cout << "  - Total steps:    " << simData.totalSteps << "\n";

    std::cout << "[PHYSICS]\n";
    std::cout << "  - Viscosity (nu): " << simData.nu << "\n\n";

    std::cout << "Initialization complete. Starting simulation... ✅\n\n";


    // Initialize steps & VTK writer
    OutputSettings output = input.output;
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
        std::cout << "Step " << simData.currStep 
                  << " done successfully in " 
                  << elapsed.count() << " s.\n";
        
        // --- ON-THE-FLY OUTPUT LOGIC ---

        // 1. Check if we should write output this step
        if (output.enabled && (simData.currStep % output.frequency == 0))
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