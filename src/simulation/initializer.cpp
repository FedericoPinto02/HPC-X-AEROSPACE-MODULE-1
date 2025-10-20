#include "simulation/initializer.hpp"
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
    const auto& grid = std::make_shared<Grid>();

    grid->nx = data.mesh.nx;
    grid->ny = data.mesh.ny;
    grid->nz = data.mesh.nz;

    grid->dx = data.mesh.dx;
    grid->dy = data.mesh.dy;
    grid->dz = data.mesh.dz;

    return grid;
}

std::shared_ptr<Field> Initializer::initializePressure(std::shared_ptr<const Grid> grid) {

    const size_t N = grid->nx * grid->ny * grid->nz;
    std::vector<double> p0(N, data.initialConditions.p0);

    const auto& pressure = std::make_shared<Field>();

    pressure->setup(grid, p0);

    return pressure;
}

std::shared_ptr<VectorField> Initializer::initializeVelocity(std::shared_ptr<const Grid> grid) {

    const size_t N = grid->nx * grid->ny * grid->nz;
    std::vector<double> u0(N, data.initialConditions.u0);
    std::vector<double> v0(N, data.initialConditions.v0);
    std::vector<double> w0(N, data.initialConditions.w0);

    const auto& velocity = std::make_shared<VectorField>();

    velocity->setup(grid, u0, v0, w0);

    return velocity;
}