#include "simulation/initializer.hpp"
#include <core/Mesh.hpp>
#include <iostream>
#include <vector>
#include <cmath>

Initializer::Initializer(const InputData& inputData) : data(inputData) {}

SimulationData Initializer::setup(){
    auto grid = buildGrid();
    auto pressure = initializePressure(grid);
    auto velocity = initializeVelocity(grid);
    auto porosity = initializePorosity(grid);

    return SimulationData{grid, pressure, velocity, porosity, data};
}

std::shared_ptr<Grid> Initializer::buildGrid() {
    auto grid = std::make_shared<Grid>(
        data.mesh.nx,
        data.mesh.ny,
        data.mesh.nz,
        data.mesh.dx,
        data.mesh.dy,
        data.mesh.dz
    );

    std::cout << "Building grid..." << std::endl;
    std::cout << "nx = " << grid->Nx << ", ny = " << grid->Ny << ", nz = " << grid->Nz << std::endl;
    std::cout << "dx = " << grid->dx << ", dy = " << grid->dy << ", dz = " << grid->dz << std::endl;

    return grid;
}

std::shared_ptr<Field> Initializer::initializePressure(std::shared_ptr<Grid> grid) {
    const size_t N = grid->Nx * grid->Ny * grid->Nz;
    std::vector<double> p0(N, data.initialConditions.p0);

    auto pressure = std::make_shared<Field>();
    pressure->setup(grid, p0);   // <-- direttamente
    return pressure;
}

std::shared_ptr<VectorField> Initializer::initializeVelocity(std::shared_ptr<Grid> grid) {
    const size_t N = grid->Nx * grid->Ny * grid->Nz;
    std::vector<double> u0(N, data.initialConditions.u0);
    std::vector<double> v0(N, data.initialConditions.v0);
    std::vector<double> w0(N, data.initialConditions.w0);

    auto velocity = std::make_shared<VectorField>();
    velocity->setup(grid, u0, v0, w0); // <-- direttamente
    return velocity;
}


std::shared_ptr<Field> Initializer::initializePorosity(std::shared_ptr<Grid> grid) {
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double dx = grid->dx;
    const double dy = grid->dy;
    const double dz = grid->dz;

    std::vector<double> kValues(Nx * Ny * Nz);

    // Calcola la porosità come funzione della posizione
    for (size_t k = 0; k < Nz; ++k) {
        for (size_t j = 0; j < Ny; ++j) {
            for (size_t i = 0; i < Nx; ++i) {
                double x = i * dx;
                double y = j * dy;
                double z = k * dz;

                // Esempio: porosità k(x,y,z) = sin(x) * sin(y) * sin(z)
                kValues[i + Nx * (j + Ny * k)] = std::sin(x) * std::sin(y) * std::sin(z);

            }
        }
    }

    auto porosity = std::make_shared<Field>();
    porosity->setup(grid, kValues);

    return porosity;
}


