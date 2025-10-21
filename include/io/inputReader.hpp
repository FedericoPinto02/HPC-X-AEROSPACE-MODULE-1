#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include <string>

/**
 * @brief Structure to hold mesh configuration data.
 */
struct MeshData {
    int nx, ny, nz;
    double dx, dy, dz;
};

/**
 * @brief Structure to hold physics parameters.
 */
struct PhysicsData {
    double viscosity;
};

/**
 * @brief Structure to hold initial conditions.
 */
struct InitialConditions {
    double u0, v0, w0, p0;
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
    InitialConditions initialConditions;
    TimeData time;
};

/**
 * @brief Class responsible for reading and parsing input configuration files.
 */
class InputReader {
public:
    /**
     * @brief Default constructor.
     */
    InputReader() = default;

    /**
     * @brief Read and parse input data from a JSON configuration file.
     * @param filename path to the configuration file
     * @return InputData structure containing all parsed data
     * @throws std::runtime_error if file cannot be opened or parsed
     */
    InputData readAndSetInput(const std::string& filename);
};

#endif // INPUTREADER_HPP