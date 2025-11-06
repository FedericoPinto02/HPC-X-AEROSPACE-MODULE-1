#include <iostream>
#include <iomanip>
#include "io/inputReader.hpp"
#include "simulation/initializer.hpp"
#include "simulation/SimulationContext.hpp"

int main() {
    try {
        // ----- 1. Lettura input -----
        InputReader reader;
        InputData input = reader.read("../data/config.json");
        std::cout << "[OK] File di input letto correttamente.\n\n";

        // ----- 2. Stampa riepilogo dei dati letti -----
        std::cout << "=== INPUT SUMMARY ===\n";

        // --- Mesh ---
        std::cout << "[MESH]\n";
        std::cout << "  Nx = " << input.mesh.nx
                  << ", Ny = " << input.mesh.ny
                  << ", Nz = " << input.mesh.nz << "\n";
        std::cout << "  dx = " << input.mesh.dx
                  << ", dy = " << input.mesh.dy
                  << ", dz = " << input.mesh.dz << "\n\n";

        // --- Time ---
        std::cout << "[TIME]\n";
        std::cout << "  dt = " << input.time.dt
                  << ", t_end = " << input.time.t_end << "\n\n";

        // --- Physics ---
        std::cout << "[PHYSICS]\n";
        std::cout << "  nu = " << input.physics.nu << "\n";
        std::cout << "  k_expr = \"" << input.physics.k_expr << "\"\n\n";

        // --- Forces ---
        std::cout << "[FORCES]\n";
        std::cout << "  fx_expr = \"" << input.forces.fx_expr << "\"\n";
        std::cout << "  fy_expr = \"" << input.forces.fy_expr << "\"\n";
        std::cout << "  fz_expr = \"" << input.forces.fz_expr << "\"\n\n";

        // --- Initial Conditions ---
        std::cout << "[INITIAL CONDITIONS]\n";
        std::cout << "  u_expr = \"" << input.initial_conditions.u_expr << "\"\n";
        std::cout << "  v_expr = \"" << input.initial_conditions.v_expr << "\"\n";
        std::cout << "  w_expr = \"" << input.initial_conditions.w_expr << "\"\n";
        std::cout << "  p_expr = \"" << input.initial_conditions.p_expr << "\"\n\n";

        // --- Boundary Conditions ---
        std::cout << "[BOUNDARY CONDITIONS]\n";
        std::cout << "  u_expr = \"" << input.boundary_conditions.u_expr << "\"\n";
        std::cout << "  v_expr = \"" << input.boundary_conditions.v_expr << "\"\n";
        std::cout << "  w_expr = \"" << input.boundary_conditions.w_expr << "\"\n\n";

        // --- Output ---
        std::cout << "[OUTPUT]\n";
        std::cout << "  enabled      = " << std::boolalpha << input.output.enabled << "\n";
        std::cout << "  dir          = \"" << input.output.dir << "\"\n";
        std::cout << "  base_filename= \"" << input.output.baseFilename << "\"\n";
        std::cout << "  frequency    = " << input.output.frequency << "\n\n";

        // --- Logging ---
        std::cout << "[LOGGING]\n";
        std::cout << "  verbose       = " << std::boolalpha << input.logging.verbose << "\n";
        std::cout << "  logToFile     = " << input.logging.logToFile << "\n";
        std::cout << "  logToConsole  = " << input.logging.logToConsole << "\n";
        std::cout << "  filename      = \"" << input.logging.filename << "\"\n";
        std::cout << "  frequency     = " << input.logging.frequency << "\n\n";

        // ----- 3. Setup simulazione -----
        Initializer init;
        SimulationData simData = init.setup(input);
        std::cout << "[OK] SimulationData inizializzato correttamente.\n\n";

        // ----- 4. Stampa riepilogo mesh simulazione -----
        std::cout << "=== MESH INFO (da SimulationData) ===\n";
        std::cout << "Grid size: " << simData.Nx << " x " << simData.Ny << " x " << simData.Nz << "\n";
        std::cout << "Spacing: dx=" << simData.dx << ", dy=" << simData.dy << ", dz=" << simData.dz << "\n";
        std::cout << "Domain length: Lx=" << simData.Lx
                  << ", Ly=" << simData.Ly
                  << ", Lz=" << simData.Lz << "\n\n";

        // ----- 5. Stampa dati temporali -----
        std::cout << "=== TIME SETTINGS ===\n";
        std::cout << "dt = " << simData.dt
                  << ", total time = " << simData.totalSimTime
                  << ", total steps = " << simData.totalSteps << "\n\n";

        // ----- 6. Stampa dati fisici -----
        std::cout << "=== PHYSICS ===\n";
        std::cout << "Viscosity (nu): " << simData.nu << "\n\n";

        std::cout << "Tutto inizializzato con successo ✅\n";

    } catch (const std::exception& e) {
        std::cerr << "❌ Errore durante l'inizializzazione: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
