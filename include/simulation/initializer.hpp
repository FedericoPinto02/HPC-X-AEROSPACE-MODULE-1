#pragma once

#include <memory>
#include <string>

#include "../data/configFunctions.hpp"
#include "io/inputReader.hpp"
#include "core/Grid.hpp"
#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"


namespace Initializer {

    /**
     * @brief Setup and initialize the SimulationData from input data.
     */
    SimulationData setup(const InputData &inputData);

    // ---------------------------------------------------------------------
    // Scalar Fields Initialization
    // ---------------------------------------------------------------------

    Field initializeFieldFromFunc(
            const double time,
            const GridPtr &grid,
            const Func &func);

    // ---------------------------------------------------------------------
    // Vector Fields Initialization
    // ---------------------------------------------------------------------

    VectorField initializeVectorFieldFromFunc(
            const double time,
            const GridPtr &grid,
            const Func &func_u,
            const Func &func_v,
            const Func &func_w);

};