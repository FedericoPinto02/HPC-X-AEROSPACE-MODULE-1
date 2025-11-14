#include "simulation/NSBSolver.hpp"

int main() {

    std::string configFile = "../data/config.json";

    NSBSolver problem(configFile);
    problem.setup();
    problem.solve();

    return 0;
}