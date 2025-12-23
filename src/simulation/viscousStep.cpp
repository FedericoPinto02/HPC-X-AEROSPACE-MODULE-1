#include "simulation/viscousStep.hpp"

ViscousStep::ViscousStep(MpiEnv &mpi, SimulationData &simData)
        : mpi(mpi), haloHandler(mpi), data_(simData) {}


void ViscousStep::setup() {
    // --- Setup velocity-like and intermediate fields -----------------------------------------------------------------
    gradP.setup(data_.grid);
    dxxEta.setup(data_.grid);
    dyyZeta.setup(data_.grid);
    dzzU.setup(data_.grid);
    xi.setup(data_.grid);

    // --- Assemble linear system scratch variables --------------------------------------------------------------------
    matrix = TridiagMat(data_.grid->Nx);
    solver_x = std::make_unique<SchurSolver>(mpi, Axis::X, matrix);
    matrix.resize(data_.grid->Ny);
    solver_y = std::make_unique<SchurSolver>(mpi, Axis::Y, matrix);
    matrix.resize(data_.grid->Nz);
    solver_z = std::make_unique<SchurSolver>(mpi, Axis::Z, matrix);

    size_t maxSysDim = std::max({data_.grid->Nx, data_.grid->Ny, data_.grid->Nz});
    rhs.reserve(maxSysDim);
    unknown_u.reserve(maxSysDim);
    unknown_v.reserve(maxSysDim);
    unknown_w.reserve(maxSysDim);
}


void ViscousStep::run() {
    computeXi();
    closeViscousStep();
}


void ViscousStep::computeXi() {
    // Prepare ingredients
    haloHandler.exchange(data_.predictor);
    haloHandler.exchange(data_.eta);
    haloHandler.exchange(data_.zeta);
    haloHandler.exchange(data_.u);
    derive.computeGradient(data_.predictor, gradP);
    derive.computeDxx(data_.eta, dxxEta);
    derive.computeDyy(data_.zeta, dyyZeta);
    derive.computeDzz(data_.u, dzzU);

    // Constant ingredients
    const double nu_val = data_.nu;
    const double dt_val = data_.dt;
    const double dt_nu_over_2_val = dt_val * nu_val * 0.5;


    // Recipie
    // g = f  - grad(p)  - nu/k * u  + nu * (dxx eta + dyy zeta + dzz u)
    // beta = 1 + dt*nu /2/k
    // xi = u + dt/beta * g

    // Let me cook
    for (Axis axis: {Axis::X, Axis::Y, Axis::Z})
    {
        size_t size = xi(axis).getData().size();
        auto f_data = data_.f(axis).getData().data();
        auto u_data = data_.u(axis).getData().data();
        auto gradP_data = gradP(axis).getData().data();
        auto dxx_data = dxxEta(axis).getData().data();
        auto dyy_data = dyyZeta(axis).getData().data();
        auto dzz_data = dzzU(axis).getData().data();
        auto inv_k_data = data_.inv_k(axis).getData().data();
        auto xi_data = xi(axis).getData().data();

        for (size_t i = 0; i < size; i++)
        {
            // Compute g
            double g_val = f_data[i]
                           - gradP_data[i]
                           - nu_val * u_data[i] * inv_k_data[i]
                           + nu_val * (dxx_data[i] + dyy_data[i] + dzz_data[i]);
            // Compute xi
            double beta_val = 1 + dt_nu_over_2_val * inv_k_data[i];
            double inv_beta_val = 1.0 / beta_val;
            xi_data[i] = u_data[i] + dt_val * inv_beta_val * g_val;
        }
    }

    size_t i, j, k;
    auto grid = data_.grid;
    auto &nx = grid->Nx;
    auto &ny = grid->Ny;
    auto &nz = grid->Nz;
    auto &time = data_.currTime;

    if (grid->hasMinBoundary(Axis::X)) {
        auto &xi_Y = xi(Axis::Y);
        auto &xi_Z = xi(Axis::Z);
        i = 0;
        for (k = 0; k < nz; k++) {
            double physical_Yz = xi_Y.getGrid().to_z(k, xi_Y.getOffset(), xi_Y.getOffsetAxis());
            double physical_Zz = xi_Z.getGrid().to_z(k, xi_Z.getOffset(), xi_Z.getOffsetAxis());    // +0.5
            for (j = 0; j < ny; j++) {
                double physical_Yy = grid->to_y(j, xi_Y.getOffset(), xi_Y.getOffsetAxis());                 // + 0.5
                xi_Y(i, j, k) = data_.bcv(0.0, physical_Yy, physical_Yz, time);

                double physical_Zy = grid->to_y(j, xi_Z.getOffset(), xi_Z.getOffsetAxis());
                xi_Z(i, j, k) = data_.bcw(0.0, physical_Zy, physical_Zz, time);
            }
        }
    }
    if (grid->hasMinBoundary(Axis::Y)) {
        auto &xi_X = xi(Axis::X);
        auto &xi_Z = xi(Axis::Z);
        j = 0;
        for (k = 0; k < nz; k++) {
            double physical_Xz = grid->to_z(k, xi_X.getOffset(), xi_X.getOffsetAxis());
            double physical_Zz = grid->to_z(k, xi_Z.getOffset(), xi_Z.getOffsetAxis());                 // +0.5
            for (i = 0; i < nx; i++) {
                double physical_Xx = xi_X.getGrid().to_x(i, xi_X.getOffset(), xi_X.getOffsetAxis());    // +0.5
                xi_X(i, j, k) = data_.bcu(physical_Xx, 0.0, physical_Xz, time);

                double physical_Zx = xi_Z.getGrid().to_x(i, xi_Z.getOffset(), xi_Z.getOffsetAxis());
                xi_Z(i, j, k) = data_.bcw(physical_Zx, 0.0, physical_Zz, time);
            }
        }
    }
    if (grid->hasMinBoundary(Axis::Z)) {
        auto &xi_X = xi(Axis::X);
        auto &xi_Y = xi(Axis::Y);
        k = 0;
        for (j = 0; j < ny; j++) {
            double physical_Xy = grid->to_y(j, xi_X.getOffset(), xi_X.getOffsetAxis());
            double physical_Yy = grid->to_y(j, xi_Y.getOffset(), xi_Y.getOffsetAxis());                 // +0.5
            for (i = 0; i < nx; i++) {
                double physical_Xx = xi_X.getGrid().to_x(i, xi_X.getOffset(), xi_X.getOffsetAxis());    // +0.5
                xi_X(i, j, k) = data_.bcu(physical_Xx, physical_Xy, 0.0, time);

                double physical_Yx = xi_Y.getGrid().to_x(i, xi_Y.getOffset(), xi_Y.getOffsetAxis());
                xi_Y(i, j, k) = data_.bcv(physical_Yx, physical_Yy, 0.0, time);
            }
        }
    }

    if (grid->hasMaxBoundary(Axis::X)) {
        auto &xi_X = xi(Axis::X);
        i = nx - 1;
        double physical_Xx = grid->to_x(i, xi_X.getOffset(), xi_X.getOffsetAxis());    // +0.5
        for (k = 0; k < nz; k++) {
            double physical_Xz = grid->to_z(k, xi_X.getOffset(), xi_X.getOffsetAxis());
            for (j = 0; j < ny; j++) {
                double physical_Xy = grid->to_y(j, xi_X.getOffset(), xi_X.getOffsetAxis());
                xi_X(i, j, k) = data_.bcu(physical_Xx, physical_Xy, physical_Xz, time);
            }
        }
    }
    if (grid->hasMaxBoundary(Axis::Y)) {
        auto &xi_Y = xi(Axis::Y);
        j = ny - 1;
        double physical_Yy = grid->to_y(j, xi_Y.getOffset(), xi_Y.getOffsetAxis());    // +0.5
        for (k = 0; k < nz; k++) {
            double physical_Yz = grid->to_z(k, xi_Y.getOffset(), xi_Y.getOffsetAxis());
            for (i = 0; i < nx; i++) {
                double physical_Yx = grid->to_x(i, xi_Y.getOffset(), xi_Y.getOffsetAxis());
                xi_Y(i, j, k) = data_.bcv(physical_Yx, physical_Yy, physical_Yz, time);
            }
        }
    }
    if (grid->hasMaxBoundary(Axis::Z)) {
        auto &xi_Z = xi(Axis::Z);
        k = nz - 1;
        double physical_Zz = grid->to_z(k, xi_Z.getOffset(), xi_Z.getOffsetAxis());    // +0.5
        for (j = 0; j < ny; j++) {
            double physical_Zy = grid->to_y(j, xi_Z.getOffset(), xi_Z.getOffsetAxis());
            for (i = 0; i < nx; i++) {
                double physical_Zx = grid->to_x(i, xi_Z.getOffset(), xi_Z.getOffsetAxis());
                xi_Z(i, j, k) = data_.bcw(physical_Zx, physical_Zy, physical_Zz, time);
            }
        }
    }
}


void ViscousStep::closeViscousStep() {
    Axis normalAxis;    // solve on the face orthogonal to normalAxis

    //==================================================================================================================
    // --- 1. Solve Eta (Direction X) ----------------------------------------------------------------------------------
    // --- when solving Eta we fill linsys with dxx derivatives
    // --- Eta.u is then solved exploiting normal Dirichlet boundary conditions
    // --- Eta.v and Eta.w are solved exploiting tangent Dirichlet boundary conditions
    // --- Differences between normal and tangent are given by staggered grid.
    //==================================================================================================================
    {
        normalAxis = Axis::X;

        size_t sysDimension = data_.grid->Nx;
        matrix.resize(sysDimension);
        rhs.resize(sysDimension);
        unknown_u.resize(sysDimension);
        unknown_v.resize(sysDimension);
        unknown_w.resize(sysDimension);

        for (long k = 0; k < data_.grid->Nz; ++k)
        {
            for (long j = 0; j < data_.grid->Ny; ++j)
            {
                assembleLocalSystem(data_.eta, xi, Axis::X, normalAxis, 0, j, k,
                                    matrix, rhs);
                solver_x->updateMatrix(matrix);
                solver_x->preprocess();
                solver_x->solve(rhs, unknown_u);

                assembleLocalSystem(data_.eta, xi, Axis::Y, normalAxis, 0, j, k,
                                    matrix, rhs);
                solver_x->updateMatrix(matrix);
                solver_x->preprocess();
                solver_x->solve(rhs, unknown_v);

                assembleLocalSystem(data_.eta, xi, Axis::Z, normalAxis, 0, j, k,
                                    matrix, rhs);
                solver_x->updateMatrix(matrix);
                solver_x->preprocess();
                solver_x->solve(rhs, unknown_w);

                double *eta_u_xLine_offset = &data_.eta(Axis::X, 0, j, k);
                double *eta_v_xLine_offset = &data_.eta(Axis::Y, 0, j, k);
                double *eta_w_xLine_offset = &data_.eta(Axis::Z, 0, j, k);
                for (size_t i = 0; i < sysDimension; ++i)
                {
                    eta_u_xLine_offset[i] = unknown_u[i];
                    eta_v_xLine_offset[i] = unknown_v[i];
                    eta_w_xLine_offset[i] = unknown_w[i];
                }
            }
        }
    }


    //==================================================================================================================
    // --- 2. Solve Zeta (Direction Y) ---------------------------------------------------------------------------------
    // --- when solving Zeta we fill linsys with dyy derivatives
    // --- Zeta.v is then solved exploiting normal Dirichlet boundary conditions
    // --- Zeta.u and Zeta.w are solved exploiting tangent Dirichlet boundary conditions
    //==================================================================================================================
    {
        normalAxis = Axis::Y;
        size_t stride = data_.zeta(Axis::X).getStride(normalAxis);

        size_t sysDimension = data_.grid->Ny;
        matrix.resize(sysDimension);
        rhs.resize(sysDimension);
        unknown_u.resize(sysDimension);
        unknown_v.resize(sysDimension);
        unknown_w.resize(sysDimension);

        for (long k = 0; k < data_.grid->Nz; ++k)
        {
            for (long i = 0; i < data_.grid->Nx; ++i)
            {
                assembleLocalSystem(data_.zeta, data_.eta, Axis::X, normalAxis, i, 0, k,
                                    matrix, rhs);
                solver_y->updateMatrix(matrix);
                solver_y->preprocess();
                solver_y->solve(rhs, unknown_u);

                assembleLocalSystem(data_.zeta, data_.eta, Axis::Y, normalAxis, i, 0, k,
                                    matrix, rhs);
                solver_y->updateMatrix(matrix);
                solver_y->preprocess();
                solver_y->solve(rhs, unknown_v);

                assembleLocalSystem(data_.zeta, data_.eta, Axis::Z, normalAxis, i, 0, k,
                                    matrix, rhs);
                solver_y->updateMatrix(matrix);
                solver_y->preprocess();
                solver_y->solve(rhs, unknown_w);

                double *zeta_u_xLine_offset = &data_.zeta(Axis::X, i, 0, k);
                double *zeta_v_xLine_offset = &data_.zeta(Axis::Y, i, 0, k);
                double *zeta_w_xLine_offset = &data_.zeta(Axis::Z, i, 0, k);
                for (size_t j = 0; j < sysDimension; ++j)
                {
                    zeta_u_xLine_offset[j * stride] = unknown_u[j];
                    zeta_v_xLine_offset[j * stride] = unknown_v[j];
                    zeta_w_xLine_offset[j * stride] = unknown_w[j];
                }
            }
        }
    }


    //==================================================================================================================
    // --- 3. Solve U (Direction Z) ------------------------------------------------------------------------------------
    // --- when solving U we fill linsys with dzz derivatives
    // --- U.w is then solved exploiting normal Dirichlet boundary conditions
    // --- U.v and U.u are solved exploiting tangent Dirichlet boundary conditions
    //==================================================================================================================
    {
        normalAxis = Axis::Z;
        size_t stride = data_.u(Axis::X).getStride(normalAxis);

        size_t sysDimension = data_.grid->Nz;
        matrix.resize(sysDimension);
        rhs.resize(sysDimension);
        unknown_u.resize(sysDimension);
        unknown_v.resize(sysDimension);
        unknown_w.resize(sysDimension);

        for (long j = 0; j < data_.grid->Ny; ++j)
        {
            for (long i = 0; i < data_.grid->Nx; ++i)
            {
                assembleLocalSystem(data_.u, data_.zeta, Axis::X, normalAxis, i, j, 0,
                                    matrix, rhs);
                solver_z->updateMatrix(matrix);
                solver_z->preprocess();
                solver_z->solve(rhs, unknown_u);

                assembleLocalSystem(data_.u, data_.zeta, Axis::Y, normalAxis, i, j, 0,
                                    matrix, rhs);
                solver_z->updateMatrix(matrix);
                solver_z->preprocess();
                solver_z->solve(rhs, unknown_v);

                assembleLocalSystem(data_.u, data_.zeta, Axis::Z, normalAxis, i, j, 0,
                                    matrix, rhs);
                solver_z->updateMatrix(matrix);
                solver_z->preprocess();
                solver_z->solve(rhs, unknown_w);

                double *u_u_xLine_offset = &data_.u(Axis::X, i, j, 0);
                double *u_v_xLine_offset = &data_.u(Axis::Y, i, j, 0);
                double *u_w_xLine_offset = &data_.u(Axis::Z, i, j, 0);
                for (size_t k = 0; k < sysDimension; ++k)
                {
                    u_u_xLine_offset[k * stride] = unknown_u[k];
                    u_v_xLine_offset[k * stride] = unknown_v[k];
                    u_w_xLine_offset[k * stride] = unknown_w[k];
                }
            }
        }
    }
}


void ViscousStep::assembleLocalSystem(
        const VectorField &eta, const VectorField &xi,
        const Axis fieldComponent, const Axis derivativeDirection,
        const size_t iStart, const size_t jStart, const size_t kStart,
        TridiagMat &matA, std::vector<double> &rhsC
) {
    if (rhsC.size() != matA.getSize()) { throw std::runtime_error("Dimension mismatch: rhsC must be size n."); }

    // --- Base line-offset of inv_k, eta and xi fields (for fast access) ----------------------------------------------
    const auto &inv_k_field = data_.inv_k(fieldComponent);
    const auto &eta_field = eta(fieldComponent);
    const auto &xi_field = xi(fieldComponent);
    const double *inv_k_line_offset = &inv_k_field(iStart, jStart, kStart);
    const double *eta_line_offset = &eta_field(iStart, jStart, kStart);
    const double *xi_line_offset = &xi_field(iStart, jStart, kStart);
    const size_t stride = eta_field.getStride(derivativeDirection);

    // --- Matrix diagonals --------------------------------------------------------------------------------------------
    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    // --- Grid spacings and coefficients ------------------------------------------------------------------------------
    const Grid &grid = eta.getGrid();
    double dx = grid.dx;
    double dy = grid.dy;
    double dz = grid.dz;

    double dCoef = 0;
    switch (derivativeDirection)
    {
        case Axis::X: dCoef = 1.0 / (dx * dx); break;
        case Axis::Y: dCoef = 1.0 / (dy * dy); break;
        case Axis::Z: dCoef = 1.0 / (dz * dz); break;
    }
    double dt_nu_over_2 = data_.dt * data_.nu * 0.5;


    //==================================================================================================================
    // --- DEFAULT SYSTEM COEFFICIENTS (same stencil; boundaries will be then overwritten) -----------------------------
    //==================================================================================================================
    for (long i = 0; i < matA.getSize(); i++)
    {
        double inv_k = inv_k_line_offset[i * stride];
        double beta = 1 + (dt_nu_over_2 * inv_k);
        double gamma = dt_nu_over_2 / beta;

        // --- MATRIX COEFFICIENTS -------------------------------------------------------------------------------------
        subdiag[i] = -gamma * dCoef;
        supdiag[i] = -gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;

        // --- RHS COEFFICIENTS ----------------------------------------------------------------------------------------
        double eta_prec_val_m1 = eta_line_offset[(i - 1) * stride];
        double eta_prec_val = eta_line_offset[i * stride];
        double eta_prec_val_p1 = eta_line_offset[(i + 1) * stride];
        double d2_eta_prec = (eta_prec_val_p1 + eta_prec_val_m1 - 2.0 * eta_prec_val) * dCoef;
        double xi_val = xi_line_offset[i * stride];
        rhsC[i] = xi_val - gamma * d2_eta_prec;
    }


    //==================================================================================================================
    // --- BOUNDARY CONDITIONS ENFORCEMENT (overwritten over default ones) ---------------------------------------------
    //==================================================================================================================
    BoundaryType boundaryType = (fieldComponent == derivativeDirection)
                                ? BoundaryType::Normal
                                : BoundaryType::Tangent;

    // --- LEFT INTERFACE ----------------------------------------------------------------------------------------------
    if (grid.hasMinBoundary(derivativeDirection))
    {
        double physical_x = grid.to_x(iStart, eta_field.getOffset(), eta_field.getOffsetAxis());
        double physical_y = grid.to_y(jStart, eta_field.getOffset(), eta_field.getOffsetAxis());
        double physical_z = grid.to_z(kStart, eta_field.getOffset(), eta_field.getOffsetAxis());

        std::function<double(double, double, double, double)> bc;
        switch (fieldComponent)
        {
            case Axis::X: bc = data_.bcu; break;
            case Axis::Y: bc = data_.bcv; break;
            case Axis::Z: bc = data_.bcw; break;
        }

        switch (boundaryType)
        {
            case BoundaryType::Normal:
            {
                // Force Wall Coordinate for Normal BC (LinearSys logic)
                if (derivativeDirection == Axis::X) physical_x = 0.0;
                if (derivativeDirection == Axis::Y) physical_y = 0.0;
                if (derivativeDirection == Axis::Z) physical_z = 0.0;

                double inv_k = inv_k_line_offset[0];
                double beta = 1.0 + (dt_nu_over_2 * inv_k);
                double gamma = dt_nu_over_2 / beta;

                // --- MATRIX BC ---------------------------------------------------------------------------------------
                diag.front() = 1.0 + 4.0 * gamma * dCoef;
                supdiag.front() = -4.0 / 3.0 * gamma * dCoef;

                // --- RHS BC ------------------------------------------------------------------------------------------
                double val_0 = xi_line_offset[0];
                double eta_0 = eta_line_offset[0];
                double eta_1 = eta_line_offset[1 * stride];
                rhsC.front() = val_0
                               + 4.0 / 3.0 * gamma * dCoef * eta_1
                               - 4.0 * gamma * dCoef * eta_0
                               + 8.0 / 3.0 * gamma * dCoef * (
                        bc(physical_x, physical_y, physical_z, data_.currTime - data_.dt) +
                        bc(physical_x, physical_y, physical_z, data_.currTime)
                );

                break;
            }
            case BoundaryType::Tangent:
            {
                diag.front() = 1.0;
                supdiag.front() = 0.0;

                rhsC.front() = bc(physical_x, physical_y, physical_z, data_.currTime);

                break;
            }
        }
    }
    else
    {
        // internal (half contribution)
        diag.front() *= 0.5;
        rhsC.front() *= 0.5;
    }

    // --- RIGHT INTERFACE ---------------------------------------------------------------------------------------------
    if (grid.hasMaxBoundary(derivativeDirection))
    {
        long endI = matA.getSize() - 1;
        double physical_x = grid.to_x(iStart + (derivativeDirection == Axis::X ? endI : 0),
                                      eta_field.getOffset(), eta_field.getOffsetAxis());
        double physical_y = grid.to_y(jStart + (derivativeDirection == Axis::Y ? endI : 0),
                                      eta_field.getOffset(), eta_field.getOffsetAxis());
        double physical_z = grid.to_z(kStart + (derivativeDirection == Axis::Z ? endI : 0),
                                      eta_field.getOffset(), eta_field.getOffsetAxis());

        std::function<double(double, double, double, double)> bc;
        switch (fieldComponent)
        {
            case Axis::X: bc = data_.bcu; break;
            case Axis::Y: bc = data_.bcv; break;
            case Axis::Z: bc = data_.bcw; break;
        }

        switch (boundaryType)
        {
            case BoundaryType::Normal:
            {
                diag.back() = 1.0;
                subdiag.back() = 0.0;

                rhsC.back() = bc(physical_x, physical_y, physical_z, data_.currTime);

                break;
            }
            case BoundaryType::Tangent:
            {
                // Enforce physical position on the actual wall
                if (derivativeDirection == Axis::X) physical_x += 0.5 * dx;
                if (derivativeDirection == Axis::Y) physical_y += 0.5 * dy;
                if (derivativeDirection == Axis::Z) physical_z += 0.5 * dz;

                double inv_k = inv_k_line_offset[endI * stride];
                double beta = 1.0 + (dt_nu_over_2 * inv_k);
                double gamma = dt_nu_over_2 / beta;

                // --- MATRIX BC ---------------------------------------------------------------------------------------
                diag.back() = 1.0 + 4.0 * gamma * dCoef;
                subdiag.back() = -4.0 / 3.0 * gamma * dCoef;

                // --- RHS BC ------------------------------------------------------------------------------------------
                double val_N = xi_line_offset[endI * stride];
                double eta_N = eta_line_offset[endI * stride];
                double eta_Nm1 = eta_line_offset[(endI - 1) * stride];
                rhsC.back() = val_N
                              + 4.0 / 3.0 * gamma * dCoef * eta_Nm1
                              - 4.0 * gamma * dCoef * eta_N
                              + 8.0 / 3.0 * gamma * dCoef * (
                                  bc(physical_x, physical_y, physical_z, data_.currTime - data_.dt) +
                                  bc(physical_x, physical_y, physical_z, data_.currTime)
                              );
            }
        }
    }
    else
    {
        diag.back() *= 0.5;
        rhsC.back() *= 0.5;
    }
}