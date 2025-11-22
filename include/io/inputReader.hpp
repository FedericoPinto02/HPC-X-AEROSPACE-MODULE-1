#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include <string>
#include <stdexcept>
#include <simulation/SimulationContext.hpp>
#include <nlohmann/json.hpp>

/**
 * @brief Structure to hold mesh configuration data.
 */
struct MeshData {
    int nx, ny, nz;
    double dx, dy, dz;
    bool input_for_manufactured_solution;
};

/**
 * @brief Structure to hold physics parameters.
 */
struct PhysicsData {
    double nu;              // kinematic viscosity
};

/**
 * @brief Structure to hold time integration parameters.
 */
struct TimeData {
    double dt;
    double t_end;
};

/**
 * @brief Structure to hold all input data from configuration file.
 */
struct InputData {
    MeshData mesh;
    PhysicsData physics;
    TimeData time;
    OutputSettings output;
    LoggingSettings logging;
    ParallelizationSettings parallelization;
};


namespace InputReader {

    /**
     * @brief Read and parse input data from a JSON configuration file.
     * @param filename path to the configuration file
     * @return InputData structure containing all parsed data
     */
    InputData read(const std::string& filename);

}

#endif // INPUTREADER_HPP
