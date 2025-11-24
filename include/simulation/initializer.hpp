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

    Field initializeFieldFromSpatialFunc(
            const GridPtr &grid,
            const Functions::Func &func);

    Field initializeFieldFromTemporalFunc(
            const double time,
            const GridPtr &grid,
            const Functions::Func &func);

    // ---------------------------------------------------------------------
    // Vector Fields Initialization
    // ---------------------------------------------------------------------

    VectorField initializeVectorFieldFromSpatialFunc(
            const GridPtr &grid,
            const Functions::Func &func_u,
            const Functions::Func &func_v,
            const Functions::Func &func_w);

    VectorField initializeVectorFieldFromTemporalFunc(
            const double time,
            const GridPtr &grid,
            const Functions::Func &func_u,
            const Functions::Func &func_v,
            const Functions::Func &func_w);

};