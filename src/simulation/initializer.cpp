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
// buildGrid
// =========================================================
std::shared_ptr<Grid> Initializer::buildGrid(const InputData& data) {
    auto grid = std::make_shared<Grid>(
        static_cast<size_t>(data.mesh.nx),
        static_cast<size_t>(data.mesh.ny),
        static_cast<size_t>(data.mesh.nz),
        data.mesh.dx,
        data.mesh.dy,
        data.mesh.dz
    );

    std::cout << "Grid built: "
              << grid->Nx << "x" << grid->Ny << "x" << grid->Nz
              << " spacing = (" << grid->dx << ", " << grid->dy << ", " << grid->dz << ")\n";
    return grid;
}

// =========================================================
// initializeFieldFromExpr
// =========================================================
Field Initializer::initializeFieldFromExpr(const std::shared_ptr<Grid>& grid, const std::string& expr) {
    auto f = makeSpatialFunc(expr);

    Field field;
    field.setup(grid, std::vector<double>(grid->size(), 0.0));
    field.populate(f);  // cell-centered di default

    return field;
}

// =========================================================
// initializeVectorFieldFromExpr
// =========================================================
VectorField Initializer::initializeVectorFieldFromExpr(
    const std::shared_ptr<Grid>& grid,
    const std::string& expr_u,
    const std::string& expr_v,
    const std::string& expr_w
) {
    auto fu = makeSpatialFunc(expr_u);
    auto fv = makeSpatialFunc(expr_v);
    auto fw = makeSpatialFunc(expr_w);

    VectorField vec;
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
    sim.Lx = sim.Nx * sim.dx;
    sim.Ly = sim.Ny * sim.dy;
    sim.Lz = sim.Nz * sim.dz;

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

    // --- Boundary conditions ---
    sim.uBoundNew = initializeVectorFieldFromExpr(grid,
                                                  inputData.boundary_conditions.u_expr,
                                                  inputData.boundary_conditions.v_expr,
                                                  inputData.boundary_conditions.w_expr);
    sim.uBoundOld = sim.uBoundNew; // copia a t=0

    // --- Permeability ---
    sim.k = initializeFieldFromExpr(grid, inputData.physics.k_expr);

    // --- Body force ---
    sim.f = initializeVectorFieldFromExpr(grid,
                                          inputData.forces.fx_expr,
                                          inputData.forces.fy_expr,
                                          inputData.forces.fz_expr);

    return sim;
}
