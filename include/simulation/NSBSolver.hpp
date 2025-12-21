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
    NSBSolver(const std::string &configFile, MpiEnv &mpi);

    void setup();
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
