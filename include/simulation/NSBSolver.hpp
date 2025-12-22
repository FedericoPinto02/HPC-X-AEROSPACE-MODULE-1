#pragma once

#include <memory>
#include <string>

#include "io/inputReader.hpp"
#include "io/VTKWriter.hpp"
#include "io/logWriter.hpp"
#include "simulation/pressureStep.hpp"
#include "simulation/viscousStep.hpp"
#include "core/Grid.hpp"
#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Class responsible for solving the Navier-Stokes-Brinkman equations.
 */
class NSBSolver
{
public:
    /**
     * @brief Constructor that initializes the solver with a configuration file and MPI environment.
     * @param configFile path to the configuration file
     * @param mpi the MPI environment
     */
    NSBSolver(const std::string &configFile, MpiEnv &mpi);

    /// Sets up the simulation by reading input data, initializing simulation data, and preparing steps and writers.
    void setup();

    /// Executes the time-stepping loop to solve the Navier-Stokes-Brinkman equations.
    void solve();

private:
    MpiEnv &mpi;

    std::string configFile;
    InputData input;
    SimulationData simData;
    OutputSettings outputSettings;
    LoggingSettings loggingSettings;

    std::unique_ptr<ViscousStep> viscousStep;
    std::unique_ptr<PressureStep> pressureStep;
    std::unique_ptr<VTKWriter> vtkWriter;
    std::unique_ptr<LogWriter> logger;
};
