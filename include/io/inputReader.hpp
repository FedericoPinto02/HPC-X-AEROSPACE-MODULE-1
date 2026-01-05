#pragma once

#include <string>
#include <stdexcept>
#include <simulation/SimulationContext.hpp>
#include <nlohmann/json.hpp>

/**
 * @struct MeshData
 * @brief Configuration parameters for the spatial discretization.
 * * This structure contains the number of nodes and the grid spacing for each 
 * Cartesian direction.
 */
struct MeshData
{
    int nx;       ///< Number of grid points in X direction.
    int ny;       ///< Number of grid points in Y direction.
    int nz;       ///< Number of grid points in Z direction.
    double dx;    ///< Grid spacing in X direction.
    double dy;    ///< Grid spacing in Y direction.
    double dz;    ///< Grid spacing in Z direction.
    /** * @brief Flag for Method of Manufactured Solutions (MMS).
     * If true, grid parameters are automatically recalculated to match MMS requirements.
     */
    bool input_for_manufactured_solution; 
};

/**
 * @struct PhysicsData
 * @brief Physical properties of the fluid.
 */
struct PhysicsData
{
    double nu;    ///< Kinematic viscosity of the fluid.
};

/**
 * @struct TimeData
 * @brief Parameters for the temporal discretization and simulation duration.
 */
struct TimeData
{
    double dt;     ///< Time step size.
    double t_end;  ///< Total simulation time.
};

/**
 * @struct InputData
 * @brief Container for all simulation configuration parameters.
 * * Aggregates mesh, physics, time, output, and logging settings parsed from 
 * the external configuration file.
 */
struct InputData
{
    MeshData mesh;          ///< Spatial discretization settings.
    PhysicsData physics;    ///< Fluid properties.
    TimeData time;          ///< Temporal integration settings.
    OutputSettings output;  ///< VTK and data export settings.
    LoggingSettings logging; ///< Console and file logging settings.
};

/**
 * @namespace InputReader
 * @brief Provides functionality to read and validate simulation parameters from JSON files.
 */
namespace InputReader
{

    /**
     * @brief Parses a JSON configuration file and populates an InputData structure.
     * * The function reads the specified file and validates the presence of mandatory sections 
     * (mesh, time, physics, output, logging).
     * @param filename Path to the `.json` configuration file.
     * @return A populated InputData structure.
     * @throws std::runtime_error If the file cannot be opened, if JSON syntax is invalid,
     * or if mandatory keys are missing.
     */
    InputData read(const std::string &filename);

}