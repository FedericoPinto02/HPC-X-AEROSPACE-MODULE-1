#include "core/MpiEnv.hpp"
#include "simulation/NSBSolver.hpp"

int main(int argc, char *argv[]) {
    std::string configFile = "../data/config.json";

    MpiEnv mpiEnv(argc, argv);

    NSBSolver problem(configFile);
    problem.setup();
    problem.solve();

    return 0;
}