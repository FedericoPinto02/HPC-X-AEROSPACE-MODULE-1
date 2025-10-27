#include <iostream>
#include <random>
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

        size_t Nx = sim.mesh->Nx, Ny = sim.mesh->Ny, Nz = sim.mesh->Nz;

        // --- Test valori porosità ---
        std::cout << "\n=== Porosity k sample values ===\n";
        std::cout << "k(0,0,0) = " << sim.porosity->operator()(5,5,5) << "\n";
        std::cout << "k(Nx/2,Ny/2,Nz/2) = " 
                << sim.porosity->operator()(Nx/2, Ny/2, Nz/2) << "\n";
        

        // --- Step 0: stampa valori iniziali e media ---
        double sum_p = 0.0, sum_u = 0.0, sum_v = 0.0, sum_w = 0.0, sum_k = 0.0;
        for (size_t i = 0; i < Nx; ++i)
            for (size_t j = 0; j < Ny; ++j)
                for (size_t k = 0; k < Nz; ++k) {
                    sum_p += sim.pressure->operator()(i,j,k);
                    sum_u += sim.velocity->x()(i,j,k);
                    sum_v += sim.velocity->y()(i,j,k);
                    sum_w += sim.velocity->z()(i,j,k);
                    sum_k += sim.porosity->operator()(i,j,k);
                }
        size_t N = Nx * Ny * Nz;
        std::cout << "\n=== Step 0: Average values ===\n";
        std::cout << "Pressure: " << sum_p / N << "\n";
        std::cout << "Velocity u: " << sum_u / N << "\n";
        std::cout << "Velocity v: " << sum_v / N << "\n";
        std::cout << "Velocity w: " << sum_w / N << "\n";
        std::cout << "Porosity k: " << sum_k / N << "\n";

        // --- Scrive il primo VTK ---
        VTKWriter writer(Nx, Ny, Nz, sim.mesh->dx, sim.mesh->dy, sim.mesh->dz);
        writer.write_legacy("../results/step_0000.vtk", sim.pressure, sim.velocity, "Initial Conditions");
        std::cout << "✅ VTK step_0000.vtk scritto.\n";

        // --- Step 1: assegna valori casuali a pressione e velocità ---
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        for (size_t i = 0; i < Nx; ++i)
            for (size_t j = 0; j < Ny; ++j)
                for (size_t k = 0; k < Nz; ++k) {
                    sim.pressure->operator()(i,j,k) = dis(gen);
                    sim.velocity->x()(i,j,k) = dis(gen);
                    sim.velocity->y()(i,j,k) = dis(gen);
                    sim.velocity->z()(i,j,k) = dis(gen);
                }

        // --- Stampa valori medi step 1 ---
        sum_p = sum_u = sum_v = sum_w = 0.0;
        for (size_t i = 0; i < Nx; ++i)
            for (size_t j = 0; j < Ny; ++j)
                for (size_t k = 0; k < Nz; ++k) {
                    sum_p += sim.pressure->operator()(i,j,k);
                    sum_u += sim.velocity->x()(i,j,k);
                    sum_v += sim.velocity->y()(i,j,k);
                    sum_w += sim.velocity->z()(i,j,k);
                }
        std::cout << "\n=== Step 1: Average values ===\n";
        std::cout << "Pressure: " << sum_p / N << "\n";
        std::cout << "Velocity u: " << sum_u / N << "\n";
        std::cout << "Velocity v: " << sum_v / N << "\n";
        std::cout << "Velocity w: " << sum_w / N << "\n";

        // --- Scrive il secondo VTK ---
        writer.write_legacy("../results/step_0001.vtk", sim.pressure, sim.velocity, "Step 1 Random Fields");
        std::cout << "✅ VTK step_0001.vtk scritto.\n";

    } catch (const std::exception& e) {
        std::cerr << "❌ Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
