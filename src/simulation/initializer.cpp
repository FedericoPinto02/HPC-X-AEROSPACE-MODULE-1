#include "simulation/initializer.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

// =========================================================
// 1. Scalar Fields Initialization
// =========================================================

Field Initializer::initializeFieldFromFunc(
        const double time,
        const GridPtr &grid,
        const Func &func
) {
    Field field;
    field.setup(grid, func);
    field.populate(time);
    return field;
}

// =========================================================
// 2. Vector Fields Initialization
// =========================================================

VectorField Initializer::initializeVectorFieldFromFunc(
        const double time,
        const GridPtr &grid,
        const Func &func_u,
        const Func &func_v,
        const Func &func_w) {
    VectorField vec;
    vec.setup(grid, func_u, func_v, func_w);
    vec.populate(time);
    return vec;
}


// =========================================================
// setup
// =========================================================
SimulationData Initializer::setup(const InputData &inputData, const MpiEnv &mpi) {
    SimulationData sim;

    // --- Grid ---
    GridPtr grid = std::make_shared<const Grid>(
            static_cast<size_t>(inputData.mesh.nx),
            static_cast<size_t>(inputData.mesh.ny),
            static_cast<size_t>(inputData.mesh.nz),
            inputData.mesh.dx, inputData.mesh.dy, inputData.mesh.dz,
            mpi);
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
    Func ui_func = ConfigFuncs::u_init_func;
    Func vi_func = ConfigFuncs::v_init_func;
    Func wi_func = ConfigFuncs::w_init_func;
    Func pi_func = ConfigFuncs::p_init_func;

    sim.u = initializeVectorFieldFromFunc(t0, grid, ui_func, vi_func, wi_func);
    sim.eta = sim.u;
    sim.zeta = sim.u;

    sim.p = initializeFieldFromFunc(t0, grid, pi_func);

    // --- Permeability ---
    Func inv_k_func = [](double x, double y, double z, double t) {
        const double INV_K_MAX_PENALIZATION = 1e10; // i.e., applied to k < 1e-10
        const double inv_k_val = 1 / ConfigFuncs::k_func(x, y, z, t);
        return inv_k_val >= INV_K_MAX_PENALIZATION ? INV_K_MAX_PENALIZATION : inv_k_val;
    };
    sim.inv_k = initializeFieldFromFunc(t0, grid, inv_k_func);

    // --- BC functions ---
    sim.bcu = ConfigFuncs::bcu_func;
    sim.bcv = ConfigFuncs::bcv_func;
    sim.bcw = ConfigFuncs::bcw_func;


    // --- Force functions ---
    sim.fx = ConfigFuncs::fx_func;
    sim.fy = ConfigFuncs::fy_func;
    sim.fz = ConfigFuncs::fz_func;

    sim.f = initializeVectorFieldFromFunc(t0, grid, sim.fx, sim.fy, sim.fz);

    return sim;
}