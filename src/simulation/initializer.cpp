#include "simulation/initializer.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>


// =========================================================
// 1. Scalar Fields Initialization
// =========================================================

Field Initializer::initializeFieldFromSpatialFunc(
    const std::shared_ptr<Grid>& grid, 
    const SpatialFunc& func
) {
    Field field;
    size_t size = grid->size();
    field.setup(grid, std::vector<double>(size, 0.0));
    field.populate(func, FieldOffset::CELL_CENTERED); 

    return field;
}

Field Initializer::initializeFieldFromTemporalFunc(
    const double time,
    const std::shared_ptr<Grid>& grid, 
    const TemporalFunc& func
) {
    Field field;
    size_t size = grid->size();
    field.setup(grid, std::vector<double>(size, 0.0));
    field.populate(time, func, FieldOffset::CELL_CENTERED); 

    return field;
}

// =========================================================
// 2. Vector Fields Initialization
// =========================================================

VectorField Initializer::initializeVectorFieldFromSpatialFunc(
    const std::shared_ptr<Grid>& grid,
    const SpatialFunc& func_u,
    const SpatialFunc& func_v,
    const SpatialFunc& func_w
) {
    VectorField vec;
    size_t size = grid->size();
    vec.setup(
        grid, 
        std::vector<double>(size, 0.0), // X component initialized to zero
        std::vector<double>(size, 0.0), // Y component initialized to zero
        std::vector<double>(size, 0.0)  // Z component initialized to zero
    );
    vec(Axis::X).populate(func_u, FieldOffset::FACE_CENTERED, Axis::X);
    vec(Axis::Y).populate(func_v, FieldOffset::FACE_CENTERED, Axis::Y);
    vec(Axis::Z).populate(func_w, FieldOffset::FACE_CENTERED, Axis::Z);

    return vec;
}

VectorField Initializer::initializeVectorFieldFromTemporalFunc(
    const double time,
    const std::shared_ptr<Grid>& grid,
    const TemporalFunc& func_u,
    const TemporalFunc& func_v,
    const TemporalFunc& func_w
) {
    VectorField vec;
    size_t size = grid->size();
    vec.setup(
        grid, 
        std::vector<double>(size, 0.0), // X component initialized to zero
        std::vector<double>(size, 0.0), // Y component initialized to zero
        std::vector<double>(size, 0.0)  // Z component initialized to zero
    );
    vec(Axis::X).populate(time, func_u, FieldOffset::FACE_CENTERED, Axis::X);
    vec(Axis::Y).populate(time, func_v, FieldOffset::FACE_CENTERED, Axis::Y);
    vec(Axis::Z).populate(time, func_w, FieldOffset::FACE_CENTERED, Axis::Z);

    return vec;
}

// =========================================================
// 3. Fields Update
// =========================================================

void Initializer::updateVectorFieldWithTemporalFunc(
    const double time,
    VectorField& vec,
    const TemporalFunc& func_u,
    const TemporalFunc& func_v,
    const TemporalFunc& func_w
) {
    vec(Axis::X).populate(time, func_u, FieldOffset::FACE_CENTERED, Axis::X);
    vec(Axis::Y).populate(time, func_v, FieldOffset::FACE_CENTERED, Axis::Y);
    vec(Axis::Z).populate(time, func_w, FieldOffset::FACE_CENTERED, Axis::Z);
}



// =========================================================
// setup
// =========================================================
SimulationData Initializer::setup(const InputData& inputData) {
    SimulationData sim;
    double t0 = 0.0;

    // --- Grid ---
    sim.Nx = static_cast<size_t>(inputData.mesh.nx);
    sim.Ny = static_cast<size_t>(inputData.mesh.ny);
    sim.Nz = static_cast<size_t>(inputData.mesh.nz);
    sim.dx = inputData.mesh.dx;
    sim.dy = inputData.mesh.dy;
    sim.dz = inputData.mesh.dz;
    sim.Lx = (sim.Nx + 0.5) * sim.dx;
    sim.Ly = (sim.Ny + 0.5) * sim.dy;
    sim.Lz = (sim.Nz + 0.5) * sim.dz;

    // --- Time ---
    sim.dt = inputData.time.dt;
    sim.currTime = t0;
    sim.totalSimTime = inputData.time.t_end;
    sim.currStep = 0;
    sim.totalSteps = static_cast<size_t>(std::ceil(sim.totalSimTime / sim.dt));

    // --- Physics ---
    sim.nu = inputData.physics.nu;

    // --- Build Grid object ---
    auto grid = std::make_shared<Grid>(sim.Nx, sim.Ny, sim.Nz, sim.dx, sim.dy, sim.dz);
    sim.gridPtr = grid; // temporary

    // --- Initial fields ---
    TemporalFunc ui_func = ConfigFuncs::u_init_func;
    TemporalFunc vi_func = ConfigFuncs::v_init_func;
    TemporalFunc wi_func = ConfigFuncs::w_init_func;
    TemporalFunc pi_func = ConfigFuncs::p_init_func;

    sim.u = initializeVectorFieldFromTemporalFunc(t0, grid, ui_func, vi_func, wi_func);
    sim.eta = sim.u;
    sim.zeta = sim.u;

    sim.p = initializeFieldFromTemporalFunc(t0, grid, pi_func);

    // --- Permeability ---
    SpatialFunc k_func = ConfigFuncs::k_func;

    sim.k = initializeFieldFromSpatialFunc(grid, k_func);

    // --- BC functions ---
    sim.bcu = ConfigFuncs::bcu_func;
    sim.bcv = ConfigFuncs::bcv_func;
    sim.bcw = ConfigFuncs::bcw_func;

    sim.uBoundNew = initializeVectorFieldFromTemporalFunc(t0, grid, sim.bcu, sim.bcv, sim.bcw);
    sim.uBoundOld = sim.uBoundNew;

    // --- Force functions ---
    sim.fx = ConfigFuncs::fx_func;
    sim.fy = ConfigFuncs::fy_func;
    sim.fz = ConfigFuncs::fz_func;

    sim.f = initializeVectorFieldFromTemporalFunc(t0, grid, sim.fx, sim.fy, sim.fz);

    return sim;
}
