#include "io/inputReader.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

InputData InputReader::read(const std::string &filename)
{
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        throw std::runtime_error("Could not open input file: " + filename);
    }

    json jsonData;
    try
    {
        inputFile >> jsonData;
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error("JSON parsing error: " + std::string(e.what()));
    }

    InputData data;

    // ----------------------------
    // Mesh
    // ----------------------------
    try
    {
        data.mesh.nx = jsonData["mesh"]["nx"];
        data.mesh.ny = jsonData["mesh"]["ny"];
        data.mesh.nz = jsonData["mesh"]["nz"];
        data.mesh.dx = jsonData["mesh"]["dx"];
        data.mesh.dy = jsonData["mesh"]["dy"];
        data.mesh.dz = jsonData["mesh"]["dz"];
        data.mesh.input_for_manufactured_solution = jsonData["mesh"]["input_for_manufactured_solution"];
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error("Error parsing 'mesh' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Time
    // ----------------------------
    try
    {
        data.time.dt = jsonData["time"]["dt"];
        data.time.t_end = jsonData["time"]["t_end"];
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error("Error parsing 'time' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Physics
    // ----------------------------
    try
    {
        data.physics.nu = jsonData["physics"]["nu"];
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error("Error parsing 'physics' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Output
    // ----------------------------
    try
    {
        data.output.enabled = jsonData["output"]["enabled"];
        data.output.dir = jsonData["output"]["dir"];
        data.output.baseFilename = jsonData["output"]["base_filename"];
        data.output.outputFrequency = jsonData["output"]["output_frequency"];
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error("Error parsing 'output' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Logging
    // ----------------------------
    try
    {
        data.logging.logToFile = jsonData["logging"]["log_to_file"];
        data.logging.logToConsole = jsonData["logging"]["log_to_console"];
        data.logging.dir = jsonData["logging"]["dir"];
        data.logging.filename = jsonData["logging"]["filename"];
    }
    catch (const json::exception &e)
    {
        throw std::runtime_error("Error parsing 'logging' section: " + std::string(e.what()));
    }

    // ----------------------------
    // Manufactured Solution Overrides
    // ----------------------------
    if (data.mesh.input_for_manufactured_solution)
    {
        // Set Characteristic Length L equal to kinematic viscosity nu.
        // This ensures that the Reynolds number Re = (U * L) / nu is exactly 1 (assuming U ~ 1).
        double L = data.physics.nu;

        // Calculate grid spacing dx based on the provided nx.
        // The domain length is L, and we use the formula: dx = L / (nx - 0.5).
        double calculated_dx = L / (static_cast<double>(data.mesh.nx) - 0.5);

        // Enforce an isotropic grid (cubic cells and cubic domain):
        // 1. Set all spatial steps equal to the calculated dx
        data.mesh.dx = calculated_dx;
        data.mesh.dy = calculated_dx;
        data.mesh.dz = calculated_dx;

        // 2. Set y and z grid points equal to nx
        data.mesh.ny = data.mesh.nx;
        data.mesh.nz = data.mesh.nx;
    }

    return data;
}