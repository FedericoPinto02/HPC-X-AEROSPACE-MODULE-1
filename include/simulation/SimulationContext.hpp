#pragma once

#include <cstddef>
#include <string>
#include <vector>
#include <functional>
#include <memory>

#include "core/Fields.hpp"

/**
 * @brief Configuration parameters for outputting simulation results (VTK).
 */
struct OutputSettings {
    std::string dir;            ///< Output directory
    std::string baseFilename;   ///< Base name for VTK files
    bool enabled;               ///< Enable/disable output
    size_t outputFrequency;     ///< Frequency of saving results
};

/**
 * @brief Configuration parameters for simulation logging.
 */
struct LoggingSettings {
    bool logToFile;             ///< Enable logging to file
    bool logToConsole;          ///< Enable logging to console
    std::string dir;            ///< Log directory
    std::string filename;       ///< Log filename
    size_t loggingFrequency;    ///< Frequency of log updates
};

/**
 * @brief Central structure containing all fields and physical properties of the simulation.
 */
struct SimulationData {
    // Grid
    GridPtr grid;               ///< Pointer to the simulation grid

    // Time integration settings
    double dt;                  ///< Time step size
    double currTime;            ///< Current simulation time
    double totalSimTime;        ///< Total time of the simulation
    size_t currStep;            ///< Current step index
    size_t totalSteps;          ///< Total number of steps

    // Velocity fields
    VectorField eta;            ///< Intermediate velocity field eta
    VectorField zeta;           ///< Intermediate velocity field zeta
    VectorField u;              ///< Final velocity field

    // Pressure field
    Field p;                    ///< Pressure field
    Field predictor;            ///< Pressure predictor field

    // Boundary conditions
    Func bcu;                   ///< BC function for u-velocity
    Func bcv;                   ///< BC function for v-velocity
    Func bcw;                   ///< BC function for w-velocity

    // Physical properties
    double nu;                  ///< Kinematic viscosity

    /** * @brief Inverse of the permeability (solid: high; fluid: low).
     * Multiplication with inv_k is preferred over division by k.
     */
    VectorField inv_k;

    // Body forces
    VectorField f;              ///< Total force field
    Func fx;                    ///< Force function in x-direction
    Func fy;                    ///< Force function in y-direction
    Func fz;                    ///< Force function in z-direction
};