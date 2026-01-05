#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "io/inputReader.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @class LogWriter
 * @brief Handles simulation output to console and log files.
 * * This class is designed for MPI parallel simulations. It ensures that only the 
 * Master process performs I/O operations to prevent stream corruption. It provides 
 * formatted output for simulation headers, real-time step progress, and final 
 * performance metrics (benchmarking).
 */
class LogWriter {
public:
    /**
     * @brief Constructs the LogWriter and initializes file streams if necessary.
     * @param isMaster Boolean flag indicating if the current MPI rank is the Master.
     * @param logSettings Configuration structure containing paths and enable flags.
     */
    LogWriter(bool isMaster, const LoggingSettings &logSettings);

    /**
     * @brief Writes the provided string to stdout and/or the log file, depending on configuration.
     * @param msg The string message to be logged.
     */
    void write(const std::string &msg);

    /**
     * @brief Renders a comprehensive overview of the simulation setup.
     * Includes details on:
     * - Simulation mode (Standard vs. Manufactured Solution)
     * - MPI Topology and domain decomposition
     * - Grid resolution (Global vs. Local) and spacing
     * - Numerical parameters (dt, end time)
     * - Output frequency settings
     * * @param mpi The MPI environment context.
     * @param input Raw input data from JSON.
     * @param simData Processed simulation state and grid information.
     * @param vtkWritten Flag indicating if the initial state was successfully exported.
     */
    void printSimulationHeader(const MpiEnv &mpi, const InputData &input, const SimulationData &simData, bool vtkWritten);

    /**
     * @brief Prints the table header for the time-stepping loop.
     * Format: [STEP | TIME | CPU TIME | STATUS]
     */
    void printStepHeader();

    /**
     * @brief Logs the progress of a specific time step.
     * @param step Current iteration index.
     * @param time Current physical time.
     * @param dt Size of the current time step.
     * @param elapsedSec Wall-clock time taken by the current step (in seconds).
     * @param isOutputStep Boolean flag indicating if a VTK file was saved in this step.
     */
    void printStepProgress(int step, double time, double dt, double elapsedSec, bool isOutputStep);

    /**
     * @brief Reports final performance statistics and simulation metrics.
     * * Calculates and displays the "Target Metric" for benchmarking:
     * \f[ \text{Metric} = \frac{\text{Total CPU Time}}{\text{Total Steps} \times \text{Total Cells}} \f]
     * * @param totalCpuTimeSec Sum of time spent in the solver loop.
     * @param meanCpuTimePerCellTimestep Normalized performance metric.
     * @param totalSteps Total number of iterations performed.
     * @param totalCells Total global number of grid cells.
     */
    void printFinalSummary(double totalCpuTimeSec, double meanCpuTimePerCellTimestep, unsigned int totalSteps,
                           const unsigned int totalCells);

private:
    bool isMaster_;        ///< True if this rank is responsible for I/O.
    bool logToFile_;       ///< Flag to enable writing to a file.
    bool logToConsole_;    ///< Flag to enable printing to stdout.
    std::string logDir_;   ///< Directory where logs are saved.
    std::string filename_; ///< Name of the log file.
    std::ofstream file_;   ///< Output file stream.

    /**
     * @brief Helper to generate visual separators (lines).
     * @param width Length of the line.
     * @param c Character to use for the line.
     * @return A string containing the separator.
     */
    std::string separator(int width = 60, char c = '=') const {
        return std::string(width, c) + "\n";
    }
};