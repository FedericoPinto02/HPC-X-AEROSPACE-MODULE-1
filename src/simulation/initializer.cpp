#include "simulation/initializer.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

// =========================================================
// 1. Scalar Fields Initialization
// =========================================================

Field Initializer::initializeFieldFromSpatialFunc(
        const GridPtr &grid,
        const Functions::Func &func
) {
    Field field;
    field.setup(grid, func);
    field.populate();
    return field;
}

Field Initializer::initializeFieldFromTemporalFunc(
        const double time,
        const GridPtr &grid,
        const Functions::Func &func
) {
    Field field;
    field.setup(grid, func);
    field.populate(time);
    return field;
}

// =========================================================
// 2. Vector Fields Initialization
// =========================================================

VectorField Initializer::initializeVectorFieldFromSpatialFunc(
        const GridPtr &grid,
        const Functions::Func &func_u,
        const Functions::Func &func_v,
        const Functions::Func &func_w
) {
    VectorField vec;
    vec.setup(grid, func_u, func_v, func_w);
    vec.populate();
    return vec;
}

VectorField Initializer::initializeVectorFieldFromTemporalFunc(
        const double time,
        const GridPtr &grid,
        const Functions::Func &func_u,
        const Functions::Func &func_v,
        const Functions::Func &func_w) {
    VectorField vec;
    vec.setup(grid, func_u, func_v, func_w);
    vec.populate(time);
    return vec;
}


// =========================================================
// setup
// =========================================================
SimulationData Initializer::setup(const InputData &inputData) {
    SimulationData sim;

    // --- Grid ---
    GridPtr grid = std::make_shared<const Grid>(
            static_cast<size_t>(inputData.mesh.nx),
            static_cast<size_t>(inputData.mesh.ny),
            static_cast<size_t>(inputData.mesh.nz),
            inputData.mesh.dx, inputData.mesh.dy, inputData.mesh.dz);
    sim.grid = grid;

    // --- Time ---
    double t0 = 0.0;
    sim.dt = inputData.time.dt;
    sim.currTime = t0;
    sim.totalSimTime = inputData.time.t_end;
    sim.currStep = 0;
    sim.totalSteps = static_cast<size_t>(std::ceil(sim.totalSimTime / sim.dt));

    // --- Physics ---
    sim.nu = inputData.physics.nu;

    // --- Initial fields ---
    Functions::Func ui_func = ConfigFuncs::u_init_func;
    Functions::Func vi_func = ConfigFuncs::v_init_func;
    Functions::Func wi_func = ConfigFuncs::w_init_func;
    Functions::Func pi_func = ConfigFuncs::p_init_func;

    sim.u = initializeVectorFieldFromTemporalFunc(t0, grid, ui_func, vi_func, wi_func);
    sim.eta = sim.u;
    sim.zeta = sim.u;

    sim.p = initializeFieldFromTemporalFunc(t0, grid, pi_func);

    // --- Permeability ---
    Functions::Func k_func = ConfigFuncs::k_func;
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