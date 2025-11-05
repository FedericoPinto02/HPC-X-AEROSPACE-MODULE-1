#include <iostream>
#include <cmath>
#include <memory>
#include <chrono>

// --- Include headers ---
#include "simulation/pressureStep.hpp" 
#include "simulation/viscousStep.hpp" 
#include "simulation/SimulationData.hpp"
#include "core/Fields.hpp"
#include "numerics/derivatives.hpp" 
#include "simulation/NSBsolver.hpp"




void NSBsolver::solve(SimulationData& simData)
{

    // ----- 1. Read input -----
    InputReader reader;
    InputData input = reader.read("../data/config.json");

    std::cout << " Input read correctly.\n";

    // ----- 2. Simulation Setup  -----
    Initializer init;
    SimulationData simData = init.setup(input);
    std::cout << "SimulationData created correctly.\n";

    // Initialize steps & VTK writer
    ViscousStep viscousStep(SimulationData simData);
    PressureStep pressureStep(SimulationData simData);
    VTKWriter vtkWriter(SimulationData simData);

    // Time integration
    for (unsigned int i = 1; i < simData.totalSteps; i++)
    {
        auto start = std::chrono::high_resolution_clock::now();  // start timer
        simData.currStep =+ 1;
        simData.currTime =+ simData.dt;
        viscousStep.run();
        pressureStep.run();
        auto end = std::chrono::high_resolution_clock::now();    // stop timer
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Step " << simData.currStep 
                  << " done successfully in " 
                  << elapsed.count() << " s.\n";
        vtkWriter.write_timestep();

    }
}