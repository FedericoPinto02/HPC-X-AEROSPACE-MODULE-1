#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "core/Fields.hpp"


struct SimulationData 
{
    // Grid
    double dx, dy, dz;
    double Lx, Ly, Lz;
    size_t Nx, Ny, Nz;

    // Time integration settings
    double dt;
    double currTime;
    double totalSimTime;
    size_t currStep;
    size_t totalSteps;

    // Velocity fields
    VectorField eta;
    VectorField zeta;
    VectorField u;

    // Pressure field
    Field p;

    // Boundary conditions
    VectorField uBoundOld;
    VectorField uBoundNew;

    // Kinematic viscosity of the fluid
    double nu;     

    // Permeability of the (porous) medium (solid: very low; fluid: very high)
    Field k;

    // Body force acting on the medium 
    VectorField f; 

    // Output & Log
    std::string outputDir;
    std::string baseFilename;
    size_t outputStepFrequency;
    bool verbose;
    bool logToFile;
    std::string logFilename;
    bool logToConsole;


};