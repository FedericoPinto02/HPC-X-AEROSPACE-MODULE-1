#include <iostream>
#include "io/inputReader.hpp"
#include "simulation/initializer.hpp"

int main(){
    InputReader inputReader;
    InputData inputData = inputReader.readAndSetInput("data/config.json");

    Initializer init(inputData);
    SimulationData sim = init.setup();

    std::cout << "vel"
                << " " << sim.velocity->x()(0,0,0) << " "
                << sim.velocity->y()(0,0,0) << " "
                << sim.velocity->z()(0,0,0) << std::endl;

    std::cout << "pres"
                << " " << (*sim.pressure)(0,0,0) << std::endl;

    return 0;
}