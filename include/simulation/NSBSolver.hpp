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
 * @class NSBSolver
 * @brief Main controller for solving the Navier-Stokes-Brinkman equations.
 * * This class manages the initialization, execution of the time-loop, 
 * and orchestration of the numerical steps (viscous and pressure), 
 * as well as logging and data output.
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
    MpiEnv &mpi;                                        ///< Reference to the MPI environment

    std::string configFile;                             ///< Path to the configuration file
    InputData input;                                    ///< Parsed input parameters
    SimulationData simData;                             ///< Global simulation state and fields
    OutputSettings outputSettings;                      ///< VTK output configuration
    LoggingSettings loggingSettings;                    ///< Logger configuration

    std::unique_ptr<ViscousStep> viscousStep;           ///< Solver for the viscous step
    std::unique_ptr<PressureStep> pressureStep;         ///< Solver for the pressure step
    std::unique_ptr<VTKWriter> vtkWriter;               ///< VTK data exporter
    std::unique_ptr<LogWriter> logger;                  ///< Simulation logger
};
