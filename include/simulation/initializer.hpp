#include <memory>
#include <io/inputReader.hpp>
#include <core/Mesh.hpp>
#include <core/Fields.hpp>

struct SimulationData{
    std::shared_ptr<Grid> mesh;
    std::shared_ptr<Field> pressure;
    std::shared_ptr<VectorField> velocity;
    InputData input;
};

class Initializer {
    public:
        Initializer(const InputData& inputData) {}
        SimulationData setup();

    private:
        InputData data;

        std::shared_ptr<Grid> buildGrid();
        std::shared_ptr<Field> initializePressure(std::shared_ptr<const Grid> grid);
        std::shared_ptr<VectorField> initializeVelocity(std::shared_ptr<const Grid> grid);
};