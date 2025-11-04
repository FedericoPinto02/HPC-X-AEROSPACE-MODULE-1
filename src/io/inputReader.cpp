#include "io/inputReader.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

InputData InputReader::read(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open input file: " + filename);
    }

    json jsonData;
    try {
        inputFile >> jsonData;
    } catch (const json::exception& e) {
        throw std::runtime_error("JSON parsing error: " + std::string(e.what()));
    }

    InputData data;

    // Parse mesh data
    data.mesh.nx = jsonData["mesh"]["nx"];
    data.mesh.ny = jsonData["mesh"]["ny"];
    data.mesh.nz = jsonData["mesh"]["nz"];

    data.mesh.dx = jsonData["mesh"]["dx"];
    data.mesh.dy = jsonData["mesh"]["dy"];
    data.mesh.dz = jsonData["mesh"]["dz"];

    // Parse physics data
    data.physics.viscosity = jsonData["physics"]["viscosity"];
    
    // Parse initial conditions
    data.initialConditions.u0 = jsonData["initial_conditions"]["u0"];
    data.initialConditions.v0 = jsonData["initial_conditions"]["v0"];
    data.initialConditions.w0 = jsonData["initial_conditions"]["w0"];
    data.initialConditions.p0 = jsonData["initial_conditions"]["p0"];

    // Parse time data
    data.time.dt = jsonData["time"]["dt"];
    data.time.t_end = jsonData["time"]["t_end"];

    return data;
}