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
};

/**
 * @brief Structure to hold physics parameters.
 */
struct PhysicsData {
    double nu;              // kinematic viscosity
    std::string k_expr;     // permeability function expression
};

/**
 * @brief Structure to hold initial conditions as string expressions.
 */
struct InitialConditions {
    std::string u_expr;
    std::string v_expr;
    std::string w_expr;
    std::string p_expr;
};

/**
 * @brief Structure to hold boundary conditions as string expressions.
 */
struct BoundaryConditions {
    std::string u_expr;
    std::string v_expr;
    std::string w_expr;
};

/**
 * @brief Structure to hold body force expressions.
 */
struct ForcesData {
    std::string fx_expr;
    std::string fy_expr;
    std::string fz_expr;
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
    InitialConditions initial_conditions;
    BoundaryConditions boundary_conditions;
    ForcesData forces;
    TimeData time;
    OutputSettings output;
    LoggingSettings logging;
    ParallelizationSettings parallelization;
};

/**
 * @brief Class responsible for reading and parsing input configuration files.
 */
class InputReader {
public:
    InputReader() = default;

    /**
     * @brief Read and parse input data from a JSON configuration file.
     * @param filename path to the configuration file
     * @return InputData structure containing all parsed data
     * @throws std::runtime_error if file cannot be opened or parsed
     */
    InputData read(const std::string& filename);
};

#endif // INPUTREADER_HPP
