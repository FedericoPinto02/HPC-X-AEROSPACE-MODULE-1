#include <iostream>
#include "io/inputReader.hpp"
#include "simulation/initializer.hpp"
#include "simulation/SimulationContext.hpp"

int main() {

    // ----- 1. Lettura input -----
    InputReader reader;
    InputData input = reader.read("../data/config.json");

    std::cout << "Input letto correttamente.\n";

    // ----- 2. Setup simulazione -----
    Initializer init;
    SimulationContext ctx = init.setup(input);
    std::cout << "SimulationContext creato correttamente.\n";

    // ----- 3. Integrazione -----
    // TimeIntegrator integrator(ctx);
    // integrator.run();

    // per accedere ai dati della simulazione:
    size_t i = 3, j = 4, k = 2;

    double u_val = ctx.state.u.x()(i,j,k);
    double v_val = ctx.state.u.y()(i,j,k);
    double w_val = ctx.state.u.z()(i,j,k);
    double p_val = ctx.state.p(i,j,k);

    std::cout << "u(" << i << "," << j << "," << k << ") = " << u_val << "\n";
    std::cout << "p(" << i << "," << j << "," << k << ") = " << p_val << "\n";

    return 0;
}
