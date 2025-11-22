#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "io/inputReader.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Handles console and file logging.
 */
class LogWriter
{
public:
    LogWriter(const LoggingSettings &logSettings);

    void write(const std::string &msg);

    void printSimulationHeader(const InputData &input, const SimulationData &simData);

    void printStepHeader();

    void printStepProgress(int step, double time, double dt, double elapsedSec, bool isOutputStep);

    void printFinalSummary(double totalCpuTimeSec, double meanCpuTimePerCellTimestep, unsigned int totalSteps,
                           const unsigned int totalCells);

private:
    bool logToFile_;
    bool logToConsole_;
    std::string logDir_;
    std::string filename_;
    std::ofstream file_;

    std::string separator(int width = 60, char c = '=') const
    {
        return std::string(width, c) + "\n";
    }
};