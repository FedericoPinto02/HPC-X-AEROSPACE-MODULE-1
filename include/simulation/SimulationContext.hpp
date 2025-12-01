#ifndef NSBSOLVER_SIMULATIONCONTEXT_HPP
#define NSBSOLVER_SIMULATIONCONTEXT_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <functional>
#include <memory>

#include "core/Fields.hpp"

/**
 * @brief Structure holding configuration parameters for outputting simulation results.
 */
struct OutputSettings {
    std::string dir;
    std::string baseFilename;
    bool enabled;
    size_t outputFrequency;
};

/**
 * @brief Structure holding configuration parameters for simulation logging (console and file).
 */
struct LoggingSettings {
    bool logToFile;
    bool logToConsole;
    std::string dir;
    std::string filename;
    size_t loggingFrequency;
};

/**
 * @brief Structure holding parameters related to parallel execution and domain decomposition.
 */
struct ParallelizationSettings {
    int schurDomains;
};

/**
 * @brief Central structure containing all transient data, fields, and physical properties of the running simulation.
 */
struct SimulationData {
    // Grid
    GridPtr grid;

    // Time integration settings
    double dt;
    double currTime;
    double totalSimTime;
    size_t currStep;
    size_t totalSteps;

    // Velocity fields
    VectorField eta;
    VectorField zeta;
    VectorField u;

    // Pressure field
    Field p;

    // Boundary conditions
    Func bcu;
    Func bcv;
    Func bcw;

    // Kinematic viscosity of the fluid
    double nu;

    // Permeability of the (porous) medium (solid: very low; fluid: very high)
    Field k;

    // Body force acting on the medium
    VectorField f;
    Func fx;
    Func fy;
    Func fz;
};

#endif //NSBSOLVER_SIMULATIONCONTEXT_HPP
