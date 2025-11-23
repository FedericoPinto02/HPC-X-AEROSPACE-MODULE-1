#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <ctime>

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

NSBSolver::NSBSolver(const std::string &configFile)
    : configFile(configFile) {}

void NSBSolver::setup()
{

    // Read input file
    input = InputReader::read(configFile);

    // Init SimulationData
    simData = Initializer::setup(input);

    // Init store settings
    outputSettings = input.output;
    loggingSettings = input.logging;
    parallelizationSettings = input.parallelization;

    // Init Steps
    viscousStep = std::make_unique<ViscousStep>(simData, parallelizationSettings);
    pressureStep = std::make_unique<PressureStep>(simData, parallelizationSettings);

    // Init Logger and Writer
    logger = std::make_unique<LogWriter>(loggingSettings);
    vtkWriter = std::make_unique<VTKWriter>(outputSettings, simData);

    // Logging
    logger->printSimulationHeader(input, simData);
}

void NSBSolver::solve()
{

    logger->printStepHeader();

    std::clock_t start_cpu_time = std::clock();

    // Time integration
    for (unsigned int i = 1; i < simData.totalSteps + 1; i++)
    {
        // 1. Update Time & Fields
        simData.uBoundOld = simData.uBoundNew;
        simData.currStep++;
        simData.currTime += simData.dt;

        simData.uBoundNew.populate(simData.currTime);
        simData.f.populate(simData.currTime);

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
            vtkWritten);
    }

    std::clock_t end_cpu_time = std::clock();
    double total_cpu_time_sec = static_cast<double>(end_cpu_time - start_cpu_time) / CLOCKS_PER_SEC;

    const unsigned int total_cells = simData.grid->size();
    double steps_times_cells = simData.totalSteps * total_cells;
    double mean_cpu_time_per_cell_timestep = total_cpu_time_sec / steps_times_cells;

    logger->printFinalSummary(
        total_cpu_time_sec,
        mean_cpu_time_per_cell_timestep,
        simData.totalSteps,
        total_cells);
}