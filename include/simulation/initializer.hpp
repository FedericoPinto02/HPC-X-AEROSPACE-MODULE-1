#pragma once

#include <memory>

#include "../data/configFunctions.hpp"
#include "core/Grid.hpp"
#include "core/Fields.hpp"
#include "io/inputReader.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @namespace Initializer
 * @brief Provides static methods to initialize the simulation environment, grids, and fields.
 * * This namespace acts as the bridge between the configuration layer and the execution layer.
 * It is responsible for:
 * - Allocating the global and local MPI grids.
 * - Calculating time-step parameters.
 * - Mapping analytical functions (initial conditions) onto the discrete grid.
 * - Setting up physical properties.
 */
namespace Initializer {

    /**
     * @brief Setup the <code>SimulationData</code> based on the provided <code>InputData</code> and MPI environment.
     * @param inputData the input data containing simulation parameters
     * @param mpi the MPI environment
     * @return the initialized simulation data
     */
    SimulationData setup(const InputData &inputData, const MpiEnv &mpi);

    /**
     * @brief Initialize a scalar <code>Field</code> from a given function.
     * @param time the time at which to evaluate the function to initialize the field
     * @param grid the grid on which to initialize the field
     * @param func the function used to initialize the field
     * @return the initialized scalar field
     */
    Field initializeFieldFromFunc(
            const double time,
            const GridPtr &grid,
            const Func &func);

    /**
     * @brief Initialize a vector <code>VectorField</code> from given functions for each component.
     * @param time the time at which to evaluate the functions to initialize the vector field
     * @param grid the grid on which to initialize the vector field
     * @param func_u the function used to initialize the x-component of the vector field
     * @param func_v the function used to initialize the y-component of the vector field
     * @param func_w the function used to initialize the z-component of the vector field
     * @return the initialized vector field
     */
    VectorField initializeVectorFieldFromFunc(
            const double time,
            const GridPtr &grid,
            const Func &func_u,
            const Func &func_v,
            const Func &func_w);
};