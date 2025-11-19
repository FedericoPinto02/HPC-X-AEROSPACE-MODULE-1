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

    // ----------------------------
    // Mesh
    // ----------------------------
    try {
        data.mesh.nx = jsonData["mesh"]["nx"];
        data.mesh.ny = jsonData["mesh"]["ny"];
        data.mesh.nz = jsonData["mesh"]["nz"];
        data.mesh.dx = jsonData["mesh"]["dx"];
        data.mesh.dy = jsonData["mesh"]["dy"];
        data.mesh.dz = jsonData["mesh"]["dz"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'mesh' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Time
    // ----------------------------
    try {
        data.time.dt = jsonData["time"]["dt"];
        data.time.t_end = jsonData["time"]["t_end"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'time' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Physics
    // ----------------------------
    try {
        data.physics.nu = jsonData["physics"]["nu"];
        data.physics.k_expr = jsonData["physics"]["k"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'physics' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Initial conditions
    // ----------------------------
    try {
        data.initial_conditions.u_expr = jsonData["initial_conditions"]["u"];
        data.initial_conditions.v_expr = jsonData["initial_conditions"]["v"];
        data.initial_conditions.w_expr = jsonData["initial_conditions"]["w"];
        data.initial_conditions.p_expr = jsonData["initial_conditions"]["p"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'initial_conditions' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Boundary conditions
    // ----------------------------
    try {
        data.boundary_conditions.u_expr = jsonData["boundary_conditions"]["u"];
        data.boundary_conditions.v_expr = jsonData["boundary_conditions"]["v"];
        data.boundary_conditions.w_expr = jsonData["boundary_conditions"]["w"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'boundary_conditions' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Forces
    // ----------------------------
    try {
        data.forces.fx_expr = jsonData["forces"]["fx"];
        data.forces.fy_expr = jsonData["forces"]["fy"];
        data.forces.fz_expr = jsonData["forces"]["fz"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'forces' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Output
    // ----------------------------
    try {
        data.output.enabled = jsonData["output"]["enabled"];
        data.output.dir = jsonData["output"]["dir"];
        data.output.baseFilename = jsonData["output"]["base_filename"];
        data.output.outputFrequency = jsonData["output"]["output_frequency"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'output' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Logging
    // ----------------------------
    try {
        data.logging.logToFile = jsonData["logging"]["log_to_file"];
        data.logging.logToConsole = jsonData["logging"]["log_to_console"];
        data.logging.dir = jsonData["logging"]["dir"];
        data.logging.filename = jsonData["logging"]["filename"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'logging' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Parallelization
    // ----------------------------
    try {
        data.parallelization.schurDomains = jsonData["parallelization"]["schur_domains"];
    } catch (const json::exception& e) {
        throw std::runtime_error("Error parsing 'parallelization' section: " + std::string(e.what()));
    }

    return data;
}
