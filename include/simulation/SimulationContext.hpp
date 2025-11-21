#ifndef NSBSOLVER_SIMULATIONCONTEXT_HPP
#define NSBSOLVER_SIMULATIONCONTEXT_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <functional>
#include <memory>

#include "core/Fields.hpp"

struct OutputSettings {
    std::string dir;
    std::string baseFilename;
    bool enabled;
    size_t outputFrequency;
};

struct LoggingSettings {
    bool logToFile;
    bool logToConsole;
    std::string dir;
    std::string filename;
    size_t loggingFrequency;
};

struct ParallelizationSettings {
    int schurDomains;
};

struct SimulationData {
    // Grid
    std::shared_ptr<const Grid> gridPtr;
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
    std::function<double (double t, double x, double y, double z)> bcu;
    std::function<double (double t, double x, double y, double z)> bcv;
    std::function<double (double t, double x, double y, double z)> bcw;

    // Kinematic viscosity of the fluid
    double nu;

    // Permeability of the (porous) medium (solid: very low; fluid: very high)
    Field k;

    // Body force acting on the medium
    VectorField f;
    std::function<double (double t, double x, double y, double z)> fx;
    std::function<double (double t, double x, double y, double z)> fy;
    std::function<double (double t, double x, double y, double z)> fz;
};

#endif //NSBSOLVER_SIMULATIONCONTEXT_HPP
