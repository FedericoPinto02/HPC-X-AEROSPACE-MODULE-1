#ifndef NSBSOLVER_INITIALIZER_HPP
#define NSBSOLVER_INITIALIZER_HPP

#include <memory>
#include <string>

#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"
#include "io/MuParserXAdapter.h"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Class responsible for initializing the simulation context.
 */
class Initializer {
public:
    Initializer() = default;

    /**
     * @brief Setup and initialize the SimulationData from input data.
     */
    SimulationData setup(const InputData& inputData);

private:
    std::shared_ptr<Grid> buildGrid(const InputData& data);

    Field initializeFieldFromExpr(const std::shared_ptr<Grid>& grid, const std::string& expr);
    VectorField initializeVectorFieldFromExpr(
        const std::shared_ptr<Grid>& grid,
        const std::string& expr_u,
        const std::string& expr_v,
        const std::string& expr_w
    );

    // helpers
    static std::function<double(double,double,double)> makeSpatialFunc(const std::string& expr);
};

#endif // NSBSOLVER_INITIALIZER_HPP
