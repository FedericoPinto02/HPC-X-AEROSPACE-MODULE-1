#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include "io/inputReader.hpp" 
#include "../simulation/SimulationContext.hpp" 

/**
 * @brief Handles console and file logging.
 */
class LogWriter {
public:
    explicit LogWriter(const LoggingSettings &logSettings);

    void printInputSummary(const InputData &input);
    void printRuntimeSummary(const SimulationData &simData);
    void printLoopProgress(const SimulationData& simData, double elapsedSec);

    ~LogWriter();

private:
    bool logToFile_;
    bool logToConsole_;
    std::ofstream file_;
    std::string filename_;

    void write(const std::string& msg);
};
