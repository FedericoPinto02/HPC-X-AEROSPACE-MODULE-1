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


void NSBSolver::solve() {

    // ----- 1. Read input -----
    InputReader reader;
    InputData input = reader.read("../data/config.json");

    std::cout << " Input read correctly.\n";

    // ----- 2. Simulation Setup  -----
    Initializer init;
    SimulationData simData = init.setup(input);
    std::cout << "SimulationData created correctly.\n";

    // Initialize steps & VTK writer
    ViscousStep viscousStep{simData};
    PressureStep pressureStep{simData};
    VTKWriter vtkWriter{simData};

    // Time integration
    for (unsigned int i = 1; i < simData.totalSteps; i++)
    {
        auto start = std::chrono::high_resolution_clock::now();  // start timer
        simData.currStep++;
        simData.currTime += simData.dt; 
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