#include "simulation/initializer.hpp"
#include <iostream>
#include <cmath>

SimulationContext Initializer::setup(const InputData& data) {
    SimulationContext ctx{
        .gridPtr = buildGrid(data),
        .timeSettings = buildTimeSettings(data),
        // .outputSettings = buildOutputSettings(data),
        // .loggingSettings = buildLoggingSettings(data),
        .constants = buildConstants(buildGrid(data), data),
        .icSettings = {},
        .bcSettings = {},
        .state = {}
    };

    return ctx;
}

std::shared_ptr<Grid> Initializer::buildGrid(const InputData& data) {
    auto grid = std::make_shared<Grid>(
        data.mesh.nx,
        data.mesh.ny,
        data.mesh.nz,
        data.mesh.dx,
        data.mesh.dy,
        data.mesh.dz
    );

    std::cout << "Grid built: "
              << grid->Nx << "x" << grid->Ny << "x" << grid->Nz
              << " spacing = (" << grid->dx << ", " << grid->dy << ", " << grid->dz << ")\n";

    return grid;
}

Field Initializer::initializePressure(const std::shared_ptr<Grid>& grid, const InputData& data) {
    const size_t N = grid->Nx * grid->Ny * grid->Nz;
    std::vector<double> values(N, data.initialConditions.p0);
    Field p;
    p.setup(grid, values);
    return p;
}

VectorField Initializer::initializeVelocity(const std::shared_ptr<Grid>& grid, const InputData& data) {
    const size_t N = grid->Nx * grid->Ny * grid->Nz;
    std::vector<double> u(N, data.initialConditions.u0);
    std::vector<double> v(N, data.initialConditions.v0);
    std::vector<double> w(N, data.initialConditions.w0);
    VectorField vel;
    vel.setup(grid, u, v, w);
    return vel;
}

Field Initializer::initializePorosity(const std::shared_ptr<Grid>& grid, const InputData& data) {
    const size_t Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;
    std::vector<double> vals(Nx * Ny * Nz);

    for (size_t k = 0; k < Nz; ++k)
        for (size_t j = 0; j < Ny; ++j)
            for (size_t i = 0; i < Nx; ++i)
                vals[i + Nx * (j + Ny * k)] =
                    std::sin(i * grid->dx) * std::sin(j * grid->dy) * std::sin(k * grid->dz);

    Field f;
    f.setup(grid, vals);
    return f;
}

Constants Initializer::buildConstants(const std::shared_ptr<Grid>& grid, const InputData& data) {
    Constants c;
    c.nu.setup(grid, std::vector<double>(grid->size(), data.physics.viscosity));
    c.rho.setup(grid, std::vector<double>(grid->size(), 1.0)); // assumiamo densità costante unitaria

    c.k = initializePorosity(grid, data);  // placeholder: could be permeability
    std::vector<double> fx(grid->size(), 0.0);
    std::vector<double> fy(grid->size(), 0.0);
    std::vector<double> fz(grid->size(), -9.81);
    c.f.setup(grid, fx, fy, fz);
    return c;
}

TimeIntegrationSettings Initializer::buildTimeSettings(const InputData& data) {
    return TimeIntegrationSettings(data.time.dt, data.time.t_end);
}

// OutputSettings Initializer::buildOutputSettings(const InputData& data) {
//     OutputSettings o;
//     o.outputDir = data.output.dir;
//     o.baseFilename = data.output.baseName;
//     o.outputStepFrequency = data.output.frequency;
//     return o;
// }

// LoggingSettings Initializer::buildLoggingSettings(const InputData& data) {
//     LoggingSettings log;
//     log.verbose = data.logging.verbose;
//     log.logToFile = data.logging.toFile;
//     log.logFilename = data.logging.filename;
//     log.logToConsole = data.logging.toConsole;
//     return log;
// }
