#include "core/MpiEnv.hpp"
#include "simulation/NSBSolver.hpp"

int main(int argc, char *argv[]) {
    std::string configFile = "../data/config.json";

    MpiEnv mpi(argc, argv);

    NSBSolver problem(configFile, mpi);
    problem.setup();
    problem.solve();

    return 0;
}