#ifndef NSBSOLVER_INITIALIZER_HPP
#define NSBSOLVER_INITIALIZER_HPP

#include <memory>
#include <string>

#include "../data/configFunctions.hpp"
#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"

using TemporalFunc = std::function<double(double, double, double, double)>;
using SpatialFunc  = std::function<double(double, double, double)>;

namespace Initializer {

    /**
     * @brief Setup and initialize the SimulationData from input data.
     */
    SimulationData setup(const InputData& inputData);

    // ---------------------------------------------------------------------
    // Scalar Fields Initialization
    // ---------------------------------------------------------------------
    
    Field initializeFieldFromSpatialFunc(
        const std::shared_ptr<Grid>& grid, 
        const SpatialFunc& func
    );

    Field initializeFieldFromTemporalFunc(
        const double time,
        const std::shared_ptr<Grid>& grid, 
        const TemporalFunc& func
    );

    // ---------------------------------------------------------------------
    // Vector Fields Initialization
    // ---------------------------------------------------------------------

    VectorField initializeVectorFieldFromSpatialFunc(
        const std::shared_ptr<Grid>& grid,
        const SpatialFunc& func_u,
        const SpatialFunc& func_v,
        const SpatialFunc& func_w
    );

    VectorField initializeVectorFieldFromTemporalFunc(
        const double time,
        const std::shared_ptr<Grid>& grid,
        const TemporalFunc& func_u,
        const TemporalFunc& func_v,
        const TemporalFunc& func_w
    );

    // ---------------------------------------------------------------------
    // Field Update
    // ---------------------------------------------------------------------

    void updateVectorFieldWithTemporalFunc(
        const double time,
        VectorField& vec,
        const TemporalFunc& func_u,
        const TemporalFunc& func_v,
        const TemporalFunc& func_w
    );

};

#endif // NSBSOLVER_INITIALIZER_HPP