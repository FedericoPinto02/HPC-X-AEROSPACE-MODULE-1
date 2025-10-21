#include "simulation/initializer.hpp"
#include <core/Mesh.hpp>
#include <iostream>
#include <vector>

Initializer::Initializer(const InputData& inputData) : data(inputData) {}

SimulationData Initializer::setup(){
    auto grid = buildGrid();
    auto pressure = initializePressure(grid);
    auto velocity = initializeVelocity(grid);

    return SimulationData{grid, pressure, velocity, data};
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
