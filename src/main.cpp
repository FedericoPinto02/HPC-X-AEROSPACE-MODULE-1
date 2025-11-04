#include <iostream>
#include "io/inputReader.hpp"
#include "simulation/initializer.hpp"
#include "simulation/SimulationContext.hpp"

int main() {

    // ----- 1. Input reading -----
    InputReader reader;
    InputData input = reader.read("../data/config.json");

    std::cout << "Input successfully read.\n";

    // ----- 2. Simulation setup -----
    Initializer init;
    SimulationContext ctx = init.setup(input);
    std::cout << "SimulationContext successfully created.\n";

    // ----- 3. Integration -----
    // TimeIntegrator integrator(ctx);
    // integrator.run();

    return 0;
}
