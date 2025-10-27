#include <iostream>
#include "io/inputReader.hpp"
#include "simulation/initializer.hpp"
#include "io/VTKWriter.hpp"

int main() {
    try {
        // 1️⃣ Legge il file di input
        InputReader reader;
        InputData data = reader.readAndSetInput(
            "../data/config.json"
        );

        // 2️⃣ Inizializza la simulazione
        Initializer init(data);
        SimulationData sim = init.setup();

        // 3️⃣ Stampa informazioni generali sulla mesh
        std::cout << "=== Mesh info ===\n";
        std::cout << "Size: " << sim.mesh->Nx << " x " 
                  << sim.mesh->Ny << " x " << sim.mesh->Nz << "\n";

        // 4️⃣ Stampa valori iniziali di pressione e velocità
        std::cout << "\n=== Initial conditions ===\n";
        std::cout << "Pressure at (0,0,0): " << sim.pressure->operator()(0,0,0) << "\n";
        std::cout << "Velocity u at (0,0,0): " << sim.velocity->x()(0,0,0) << "\n";
        std::cout << "Velocity v at (0,0,0): " << sim.velocity->y()(0,0,0) << "\n";
        std::cout << "Velocity w at (0,0,0): " << sim.velocity->z()(0,0,0) << "\n";

        // 5️⃣ Calcolo dei valori medi su tutta la mesh
        double sum_p = 0.0, sum_u = 0.0, sum_v = 0.0, sum_w = 0.0;
        size_t Nx = sim.mesh->Nx, Ny = sim.mesh->Ny, Nz = sim.mesh->Nz;

        for (size_t i = 0; i < Nx; ++i)
            for (size_t j = 0; j < Ny; ++j)
                for (size_t k = 0; k < Nz; ++k) {
                    sum_p += sim.pressure->operator()(i,j,k);
                    sum_u += sim.velocity->x()(i,j,k);
                    sum_v += sim.velocity->y()(i,j,k);
                    sum_w += sim.velocity->z()(i,j,k);
                }

        size_t N = Nx * Ny * Nz;
        std::cout << "\n=== Average values over the mesh ===\n";
        std::cout << "Average pressure: " << sum_p / N << "\n";
        std::cout << "Average velocity u: " << sum_u / N << "\n";
        std::cout << "Average velocity v: " << sum_v / N << "\n";
        std::cout << "Average velocity w: " << sum_w / N << "\n";

        // 6️⃣ Scrive un file VTK con i campi iniziali
        VTKWriter writer(Nx, Ny, Nz, sim.mesh->dx, sim.mesh->dy, sim.mesh->dz);
        writer.write_legacy("../results/initial.vtk", sim.pressure, sim.velocity, "Initial Conditions");

        std::cout << "\n✅ VTK file 'results/initial.vtk' scritto correttamente!\n";

    } catch (const std::exception& e) {
        std::cerr << "❌ Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
