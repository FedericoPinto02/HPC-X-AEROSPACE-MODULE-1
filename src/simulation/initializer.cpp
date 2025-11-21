#include "simulation/initializer.hpp"

#include <iostream>
#include <cmath>
#include <algorithm> // std::fill

using namespace MuParserXAdapter; // per createFunction / createTimeFunction

// =========================================================
// makeSpatialFunc
// =========================================================
std::function<double(double,double,double)> Initializer::makeSpatialFunc(const std::string& expr) {
    if (expr.empty()) {
        return [](double, double, double) { return 0.0; };
    }

    // Se la stringa contiene 't', crea una funzione temporale e valuta a t=0
    if (expr.find('t') != std::string::npos) {
        auto ft = createTimeFunction(expr);
        return [ft](double x, double y, double z) {
            return ft(0.0, x, y, z);
        };
    } else {
        return createFunction(expr);
    }
}

// =========================================================
// makeTemporalFunc (Nuovo)
// =========================================================
TemporalFunc Initializer::makeTemporalFunc(const std::string& expr) {
    if (expr.empty()) {
        // Se l'espressione è vuota, la funzione restituisce 0.0 per ogni t, x, y, z
        return [](double, double, double, double) { return 0.0; };
    }
    // MuParserXAdapter::createTimeFunction crea una TemporalFunc (t, x, y, z)
    return createTimeFunction(expr);
}// =========================================================
// initializeFieldFromExpr (Corretto la Grid)
// =========================================================
Field Initializer::initializeFieldFromExpr(const std::shared_ptr<const Grid>& grid, const std::string& expr) {
    auto f = makeSpatialFunc(expr);

    Field field;
    // La dimensione è 'grid->size()', che è il numero di celle (CELL_CENTERED)
    field.setup(grid, std::vector<double>(grid->size(), 0.0));
    field.populate(f, FieldOffset::CELL_CENTERED);

    return field;
}

// =========================================================
// initializeVectorFieldFromExpr (Corretto la Grid)
// =========================================================
VectorField Initializer::initializeVectorFieldFromExpr(
    const std::shared_ptr<const Grid>& grid,
    const std::string& expr_u,
    const std::string& expr_v,
    const std::string& expr_w
) {
    auto fu = makeSpatialFunc(expr_u);
    auto fv = makeSpatialFunc(expr_v);
    auto fw = makeSpatialFunc(expr_w);

    VectorField vec;
    // La dimensione è 'grid->size()', che è il numero di celle (CELL_CENTERED)
    vec.setup(grid,
              std::vector<double>(grid->size(), 0.0),
              std::vector<double>(grid->size(), 0.0),
              std::vector<double>(grid->size(), 0.0));

    // populate ogni componente, usando FACE_CENTERED per i vettori
    vec(Axis::X).populate(fu, FieldOffset::FACE_CENTERED, Axis::X);
    vec(Axis::Y).populate(fv, FieldOffset::FACE_CENTERED, Axis::Y);
    vec(Axis::Z).populate(fw, FieldOffset::FACE_CENTERED, Axis::Z);

    return vec;
}

// =========================================================
// initializeVectorFieldFromTemporalFunc (Implementazione)
// =========================================================
VectorField Initializer::initializeVectorFieldFromTemporalFunc(
    const double time,
    const std::shared_ptr<const Grid>& grid,
    const TemporalFunc& func_u,
    const TemporalFunc& func_v,
    const TemporalFunc& func_w
) {
    VectorField vec;
    vec.setup(grid,
              std::vector<double>(grid->size(), 0.0),
              std::vector<double>(grid->size(), 0.0),
              std::vector<double>(grid->size(), 0.0));

    // Si assume che VectorField::Field::populate abbia un overload per funzioni temporali.
    // Basandomi sul tuo codice, si assume che sia:
    // field.populate(time, temporal_func, offset, axis)

    // Inizializza ogni componente al tempo specificato
    vec(Axis::X).populate(time, func_u, FieldOffset::FACE_CENTERED, Axis::X);
    vec(Axis::Y).populate(time, func_v, FieldOffset::FACE_CENTERED, Axis::Y);
    vec(Axis::Z).populate(time, func_w, FieldOffset::FACE_CENTERED, Axis::Z);

    return vec;
}

// =========================================================
// updateVectorFieldWithTemporalFunc (Implementazione)
// =========================================================
void Initializer::updateVectorFieldWithTemporalFunc(
    const double time,
    VectorField& vec,
    const TemporalFunc& func_u,
    const TemporalFunc& func_v,
    const TemporalFunc& func_w
) {
    // Si assume che Field::populate abbia un overload che accetta il tempo
    // e aggiorna il campo in-place.
    // field.populate(time, temporal_func, offset, axis)
    
    vec(Axis::X).populate(time, func_u, FieldOffset::FACE_CENTERED, Axis::X);
    vec(Axis::Y).populate(time, func_v, FieldOffset::FACE_CENTERED, Axis::Y);
    vec(Axis::Z).populate(time, func_w, FieldOffset::FACE_CENTERED, Axis::Z);
}



// =========================================================
// setup
// =========================================================
SimulationData Initializer::setup(const InputData& inputData) {
    SimulationData sim;

    // --- Grid / mesh ---
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
    sim.currTime = 0.0;
    sim.totalSimTime = inputData.time.t_end;
    sim.currStep = 0;
    sim.totalSteps = static_cast<size_t>(std::ceil(sim.totalSimTime / sim.dt));

    // --- Physics ---
    sim.nu = inputData.physics.nu;

    // --- Build Grid object ---
    auto grid = std::make_shared<Grid>(sim.Nx, sim.Ny, sim.Nz, sim.dx, sim.dy, sim.dz);
    sim.gridPtr = grid; // temporary

    // --- Initial fields ---
    sim.u = initializeVectorFieldFromExpr(grid,
                                          inputData.initial_conditions.u_expr,
                                          inputData.initial_conditions.v_expr,
                                          inputData.initial_conditions.w_expr);
    sim.eta = sim.u;
    sim.zeta = sim.u;

    sim.p = initializeFieldFromExpr(grid, inputData.initial_conditions.p_expr);

    // --- Permeability ---
    sim.k = initializeFieldFromExpr(grid, inputData.physics.k_expr);

    // Salva le funzioni temporali di Boundary Condition in simData
    sim.bcu = makeTemporalFunc(inputData.boundary_conditions.u_expr);
    sim.bcv = makeTemporalFunc(inputData.boundary_conditions.v_expr);
    sim.bcw = makeTemporalFunc(inputData.boundary_conditions.w_expr);

    // Salva le funzioni temporali delle Forze in simData
    sim.fx = makeTemporalFunc(inputData.forces.fx_expr);
    sim.fy = makeTemporalFunc(inputData.forces.fy_expr);
    sim.fz = makeTemporalFunc(inputData.forces.fz_expr);

    // Boundary conditions - Inizializza uBoundNew a t=0 usando le funzioni temporali
    sim.uBoundNew = initializeVectorFieldFromTemporalFunc(sim.currTime, grid,
                                                          sim.bcu, sim.bcv, sim.bcw);
    sim.uBoundOld = sim.uBoundNew; // copia a t=0

    // Body force - Inizializza f a t=0 usando le funzioni temporali
    sim.f = initializeVectorFieldFromTemporalFunc(sim.currTime, grid,
                                                  sim.fx, sim.fy, sim.fz);

    return sim;
}
