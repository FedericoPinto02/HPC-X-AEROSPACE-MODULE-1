#ifndef NSBSOLVER_INITIALIZER_HPP
#define NSBSOLVER_INITIALIZER_HPP

#include <memory>
#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Class responsible for initializing the simulation context.
 */
class Initializer {
public:
    Initializer() = default;

    /**
     * @brief Setup and initialize the simulation context from input data.
     * @param inputData configuration data read from JSON
     * @return fully constructed SimulationContext
     */
    SimulationContext setup(const InputData& inputData);

private:
    /**
     * @brief Build the computational grid from input data.
     */
    std::shared_ptr<Grid> buildGrid(const InputData& data);

    /**
     * @brief Initialize the physical fields.
     */
    Field initializePressure(const std::shared_ptr<Grid>& grid, const InputData& data);
    VectorField initializeVelocity(const std::shared_ptr<Grid>& grid, const InputData& data);
    Field initializePorosity(const std::shared_ptr<Grid>& grid, const InputData& data);

    /**
     * @brief Build constant fields (nu, rho, k, f) from input data.
     */
    Constants buildConstants(const std::shared_ptr<Grid>& grid, const InputData& data);

    /**
     * @brief Build time integration settings.
     */
    TimeIntegrationSettings buildTimeSettings(const InputData& data);

    /**
     * @brief Build output/logging settings.
     */
    // OutputSettings buildOutputSettings(const InputData& data);
    // LoggingSettings buildLoggingSettings(const InputData& data);
};

#endif // NSBSOLVER_INITIALIZER_HPP
