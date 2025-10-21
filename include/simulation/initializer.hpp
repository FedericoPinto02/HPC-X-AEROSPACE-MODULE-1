#ifndef NSBSOLVER_INITIALIZER_HPP
#define NSBSOLVER_INITIALIZER_HPP

#include <memory>
#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"

/**
 * @brief Structure to hold all simulation data.
 */
struct SimulationData {
    std::shared_ptr<Grid> mesh;
    std::shared_ptr<Field> pressure;
    std::shared_ptr<VectorField> velocity;
    InputData input;
};

/**
 * @brief Class responsible for initializing simulation data.
 */
class Initializer {
public:
    /**
     * @brief Constructor.
     * @param inputData the input data from configuration file
     */
    Initializer(const InputData& inputData);

    /**
     * @brief Setup and initialize all simulation components.
     * @return SimulationData structure containing all initialized data
     */
    SimulationData setup();

private:
    InputData data;

    /**
     * @brief Build the computational grid from input data.
     * @return shared pointer to the Grid
     */
    std::shared_ptr<Grid> buildGrid();

    /**
     * @brief Initialize the pressure field.
     * @param grid the computational grid
     * @return shared pointer to the pressure Field
     */
    std::shared_ptr<Field> initializePressure(std::shared_ptr<Grid> grid);

    /**
     * @brief Initialize the velocity vector field.
     * @param grid the computational grid
     * @return shared pointer to the velocity VectorField
     */
    std::shared_ptr<VectorField> initializeVelocity(std::shared_ptr<Grid> grid);
};

#endif // NSBSOLVER_INITIALIZER_HPP