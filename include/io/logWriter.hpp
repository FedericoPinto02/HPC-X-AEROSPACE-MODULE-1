#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "io/inputReader.hpp"
#include "simulation/SimulationContext.hpp"

/// Class for writing log messages to console and/or file.
class LogWriter {
public:
    /**
     * @brief Constructs the <code>LogWriter</code> with specified logging settings.
     * @param logSettings the logging settings
     */
    LogWriter(bool isMaster, const LoggingSettings &logSettings);

    /**
     * @brief Writes a message to the log (console and/or file).
     * @param msg the message to write
     */
    void write(const std::string &msg);

    /**
     * @brief Prints the simulation header with configuration details.
     * @param input the input data
     * @param simData the simulation data
     * @param vtkWritten whether the first VTK output has been written
     */
    void printSimulationHeader(const MpiEnv &mpi, const InputData &input, const SimulationData &simData, bool vtkWritten);

    /// Prints the header for the step progress table.
    void printStepHeader();

    /**
     * @brief Prints the progress of a simulation step.
     * @param step the current step number
     * @param time the current simulation time
     * @param dt the time step size
     * @param elapsedSec the elapsed CPU time in seconds for this step
     * @param isOutputStep whether this step included output writing
     */
    void printStepProgress(int step, double time, double dt, double elapsedSec, bool isOutputStep);

    /**
     * @brief Prints the final summary of the simulation.
     * @param totalCpuTimeSec the total CPU time in seconds
     * @param meanCpuTimePerCellTimestep the mean CPU time per cell per timestep
     * @param totalSteps the total number of steps executed
     * @param totalCells the total number of cells in the simulation
     */
    void printFinalSummary(double totalCpuTimeSec, double meanCpuTimePerCellTimestep, unsigned int totalSteps,
                           const unsigned int totalCells);

private:
    bool isMaster_;

    bool logToFile_;
    bool logToConsole_;
    std::string logDir_;
    std::string filename_;
    std::ofstream file_;

    std::string separator(int width = 60, char c = '=') const {
        return std::string(width, c) + "\n";
    }
};