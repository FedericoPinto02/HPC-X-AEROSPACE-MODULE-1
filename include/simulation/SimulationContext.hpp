#ifndef SIMULATIONCONTEXT_HPP
#define SIMULATIONCONTEXT_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "core/Fields.hpp"


struct TimeIntegrationSettings {
    double dt;
    double currTime;
    double totalTime;
    size_t currStep;
    size_t totalSteps;
};

struct OutputSettings {
    std::string outputDir;
    std::string baseFilename;
    // TODO - etc... TBD
    size_t outputStepFrequency;
};

struct LoggingSettings {
    bool verbose;
    bool logToFile;
    std::string logFilename;
    bool logToConsole;
    // TODO - etc... TBD
};

struct Constants {
    // Kinematic viscosity of the fluid (solid: very high; fluid: low)
    Field<> nu;

    // Density of the medium (solid: very high; fluid: low)
    Field<> rho;

    // Permeability of the (porous) medium (solid: very low; fluid: very high)
    Field<> k;

    // Body force acting on the medium (e.g. gravity)
    VectorField<> f;
};

struct InitialConditionsSettings {
    // TODO - TBD
};

struct BoundaryConditionsSettings {
    // TODO - TBD
};

struct SimulationState {
    /*
     * TODO - WORK IN PROGRESS... (VectorFieldHandler, FieldHandler)
    // Velocity fields
    VectorFieldHandler<> xi;
    VectorFieldHandler<> eta;
    VectorFieldHandler<> zeta;
    VectorFieldHandler<> u;

    // Pressure fields
    FieldHandler<> psi;
    FieldHandler<> phi;
    FieldHandler<> c_phi;
    FieldHandler<> p;
     */
};

struct SimulationContext {
    std::shared_ptr<const Grid> gridPtr; // todo - move into proper Mesh object ??
    TimeIntegrationSettings timeSettings;
    OutputSettings outputSettings;
    LoggingSettings loggingSettings;
    Constants constants;
    InitialConditionsSettings icSettings;
    BoundaryConditionsSettings bcSettings;
    SimulationState state; // mutable state of the simulation
};

#endif // SIMULATIONCONTEXT_HPP