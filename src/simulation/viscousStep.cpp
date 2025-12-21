#include "simulation/viscousStep.hpp"

ViscousStep::ViscousStep(MpiEnv &mpi, SimulationData &simData)
        : mpi(mpi), data_(simData) {
    g.setup(data_.grid);
    gradP.setup(data_.grid);
    dxxEta.setup(data_.grid);
    dyyZeta.setup(data_.grid);
    dzzU.setup(data_.grid);
    xi.setup(data_.grid);
}


void ViscousStep::run() {
    HaloHandler haloHandler(mpi);
    haloHandler.exchange(data_.predictor);
    haloHandler.exchange(data_.eta);
    haloHandler.exchange(data_.zeta);
    haloHandler.exchange(data_.u);
    computeG();
    computeXi();
    closeViscousStep();
}


void ViscousStep::computeG() {
    // Ingredients list
    auto &eta = data_.eta;
    auto &zeta = data_.zeta;
    auto &u = data_.u;
    auto &predictor = data_.predictor;
    double nu_val = data_.nu;

    Derivatives derive;
    derive.computeGradient(predictor, gradP);
    derive.computeDxx(eta, dxxEta);
    derive.computeDyy(zeta, dyyZeta);
    derive.computeDzz(u, dzzU);

    // Recipie
    // g = f  - grad(p)  - nu/k * u  + nu * (dxx eta + dyy zeta + dzz u)

    // Let me cook
    for (Axis axis: {Axis::X, Axis::Y, Axis::Z}) {
        auto &f_data = data_.f(axis).getData();
        auto &u_data = data_.u(axis).getData();
        auto &gradP_data = gradP(axis).getData();
        auto &dxx_data = dxxEta(axis).getData();
        auto &dyy_data = dyyZeta(axis).getData();
        auto &dzz_data = dzzU(axis).getData();
        auto &g_data = g(axis).getData();
        auto &inv_k_data = data_.inv_k(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++) {
            g_data[i] = f_data[i]
                        - gradP_data[i]
                        - nu_val * u_data[i] * inv_k_data[i]
                        + nu_val * (dxx_data[i] + dyy_data[i] + dzz_data[i]);
        }
    }
}

void ViscousStep::computeXi() {
    // Ingredients list
    const double nu_val = data_.nu;
    const double dt_val = data_.dt;
    const double dt_nu_over_2_val = dt_val * nu_val * 0.5;

    // Recipie
    // beta = 1+ dt*nu /2/k
    // xi = u + dt/beta * g

    // Let me cook
    for (Axis axis: {Axis::X, Axis::Y, Axis::Z}) {

        auto &u_data = data_.u(axis).getData();
        auto &inv_k_data = data_.inv_k(axis).getData();
        auto &g_data = g(axis).getData();
        auto &xi_data = xi(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++) {
            double beta_val = 1 + dt_nu_over_2_val * inv_k_data[i];
            xi_data[i] = u_data[i]
                         + dt_val * g_data[i]
                           / beta_val;
        }
    }

    size_t i, j, k;
    auto grid = data_.grid;
    auto &nx = grid->Nx;
    auto &ny = grid->Ny;
    auto &nz = grid->Nz;
    auto &dx = grid->dx;
    auto &dy = grid->dy;
    auto &dz = grid->dz;
    auto &time = data_.currTime;

    if (grid->hasMinBoundary(Axis::X)) {
        auto &xi_Y = xi(Axis::Y);
        auto &xi_Z = xi(Axis::Z);
        i = 0;
        for (j = 0; j < ny; j++) {
            double physical_Yy = grid->to_y(j, xi_Y.getOffset(), xi_Y.getOffsetAxis());                 // + 0.5
            double physical_Zy = grid->to_y(j, xi_Z.getOffset(), xi_Z.getOffsetAxis());
            for (k = 0; k < nz; k++) {
                double physical_Yz = xi_Y.getGrid().to_z(k, xi_Y.getOffset(), xi_Y.getOffsetAxis());
                xi_Y(i, j, k) = data_.bcv(0.0, physical_Yy, physical_Yz, time);

                double physical_Zz = xi_Z.getGrid().to_z(k, xi_Z.getOffset(), xi_Z.getOffsetAxis());    // +0.5
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


    // ------------------------------------------
    // Solve Eta --------------------------------
    // ------------------------------------------
    {
        // when solving Eta we fill linsys with dxx derivatives
        // Eta.u is then solved exploiting normal Dirichlet boundary conditions
        // Eta.v and Eta.w are solved exploiting tangent Dirichlet boundary conditions
        // Differences between normal and tangent are given by staggered grid.

        normalAxis = Axis::X;
        size_t sysDimension = data_.grid->Nx; // dimension of linear system to solve

        TridiagMat matrix_u(sysDimension);
        TridiagMat matrix_v(sysDimension);
        TridiagMat matrix_w(sysDimension);
        std::vector<double> rhs_u(sysDimension);
        std::vector<double> rhs_v(sysDimension);
        std::vector<double> rhs_w(sysDimension);
        std::vector<double> unknown_u(sysDimension);
        std::vector<double> unknown_v(sysDimension);
        std::vector<double> unknown_w(sysDimension);

        for (size_t k = 0; k < data_.grid->Nz; ++k)
        {
            for (size_t j = 0; j < data_.grid->Ny; ++j)
            {
                assembleLocalSystem(data_, data_.eta, xi, Axis::X, Axis::X, 0, j, k,
                                    matrix_u, rhs_u);
                SchurSolver solver_u(mpi, normalAxis, matrix_u);
                solver_u.preprocess();
                solver_u.solve(rhs_u, unknown_u);

                assembleLocalSystem(data_, data_.eta, xi, Axis::Y, Axis::X, 0, j, k,
                                    matrix_v, rhs_v);
                SchurSolver solver_v(mpi, normalAxis, matrix_v);
                solver_v.preprocess();
                solver_v.solve(rhs_v, unknown_v);

                assembleLocalSystem(data_, data_.eta, xi, Axis::Z, Axis::X, 0, j, k,
                                    matrix_w, rhs_w);
                SchurSolver solver_w(mpi, normalAxis, matrix_w);
                solver_w.preprocess();
                solver_w.solve(rhs_w, unknown_w);

                for (size_t i = 0; i < sysDimension; ++i)
                {
                    data_.eta(Axis::X, i, j, k) = unknown_u[i];
                    data_.eta(Axis::Y, i, j, k) = unknown_v[i];
                    data_.eta(Axis::Z, i, j, k) = unknown_w[i];
                }
            }
        }
    }



    // ------------------------------------------
    // Solve Zeta --------------------------------
    // ------------------------------------------
    {
        // when solving Zeta we fill linsys with dyy derivatives
        // Zeta.v is then solved exploiting normal Dirichlet boundary conditions
        // Zeta.u and Zeta.w are solved exploiting tangent Dirichlet boundary conditions

        normalAxis = Axis::Y;
        size_t sysDimension = data_.grid->Ny; // dimension of linear system to solve

        TridiagMat matrix_u(sysDimension);
        TridiagMat matrix_v(sysDimension);
        TridiagMat matrix_w(sysDimension);
        std::vector<double> rhs_u(sysDimension);
        std::vector<double> rhs_v(sysDimension);
        std::vector<double> rhs_w(sysDimension);
        std::vector<double> unknown_u(sysDimension);
        std::vector<double> unknown_v(sysDimension);
        std::vector<double> unknown_w(sysDimension);

        for (size_t k = 0; k < data_.grid->Nz; ++k)
        {
            for (size_t i = 0; i < data_.grid->Nx; ++i)
            {
                assembleLocalSystem(data_, data_.zeta, data_.eta, Axis::X, Axis::Y, i, 0, k,
                                    matrix_u, rhs_u);
                SchurSolver solver_u(mpi, normalAxis, matrix_u);
                solver_u.preprocess();
                solver_u.solve(rhs_u, unknown_u);

                assembleLocalSystem(data_, data_.zeta, data_.eta, Axis::Y, Axis::Y, i, 0, k,
                                    matrix_v, rhs_v);
                SchurSolver solver_v(mpi, normalAxis, matrix_v);
                solver_v.preprocess();
                solver_v.solve(rhs_v, unknown_v);

                assembleLocalSystem(data_, data_.zeta, data_.eta, Axis::Z, Axis::Y, i, 0, k,
                                    matrix_w, rhs_w);
                SchurSolver solver_w(mpi, normalAxis, matrix_w);
                solver_w.preprocess();
                solver_w.solve(rhs_w, unknown_w);

                for (size_t j = 0; j < sysDimension; ++j)
                {
                    data_.zeta(Axis::X, i, j, k) = unknown_u[j];
                    data_.zeta(Axis::Y, i, j, k) = unknown_v[j];
                    data_.zeta(Axis::Z, i, j, k) = unknown_w[j];
                }
            }
        }
    }



    // ------------------------------------------
    // Solve U --------------------------------
    // ------------------------------------------
    {
        // when solving U we fill linsys with dzz derivatives
        // U.w is then solved exploiting normal Dirichlet boundary conditions
        // U.v and U.u are solved exploiting tangent Dirichlet boundary conditions

        normalAxis = Axis::Z;
        size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve

        TridiagMat matrix_u(sysDimension);
        TridiagMat matrix_v(sysDimension);
        TridiagMat matrix_w(sysDimension);
        std::vector<double> rhs_u(sysDimension);
        std::vector<double> rhs_v(sysDimension);
        std::vector<double> rhs_w(sysDimension);
        std::vector<double> unknown_u(sysDimension);
        std::vector<double> unknown_v(sysDimension);
        std::vector<double> unknown_w(sysDimension);

        for (size_t j = 0; j < data_.grid->Ny; ++j)
        {
            for (size_t i = 0; i < data_.grid->Nx; ++i)
            {
                assembleLocalSystem(data_, data_.u, data_.zeta, Axis::X, Axis::Z, i, j, 0,
                                    matrix_u, rhs_u);
                SchurSolver solver_u(mpi, normalAxis, matrix_u);
                solver_u.preprocess();
                solver_u.solve(rhs_u, unknown_u);

                assembleLocalSystem(data_, data_.u, data_.zeta, Axis::Y, Axis::Z, i, j, 0,
                                    matrix_v, rhs_v);
                SchurSolver solver_v(mpi, normalAxis, matrix_v);
                solver_v.preprocess();
                solver_v.solve(rhs_v, unknown_v);

                assembleLocalSystem(data_, data_.u, data_.zeta, Axis::Z, Axis::Z, i, j, 0,
                                    matrix_w, rhs_w);
                SchurSolver solver_w(mpi, normalAxis, matrix_w);
                solver_w.preprocess();
                solver_w.solve(rhs_w, unknown_w);

                for (size_t k = 0; k < sysDimension; ++k)
                {
                    data_.u(Axis::X, i, j, k) = unknown_u[k];
                    data_.u(Axis::Y, i, j, k) = unknown_v[k];
                    data_.u(Axis::Z, i, j, k) = unknown_w[k];
                }
            }
        }
    }
}


void ViscousStep::assembleLocalSystem(
        const SimulationData &simData, const VectorField &eta, const VectorField &xi,
        const Axis fieldComponent, const Axis derivativeDirection,
        const size_t iStart, const size_t jStart, const size_t kStart,
        TridiagMat &matA, std::vector<double> &rhsC
) {

    std::vector<double> &diag = matA.getDiag(0);
    std::vector<double> &subdiag = matA.getDiag(-1);
    std::vector<double> &supdiag = matA.getDiag(1);

    if (rhsC.size() != matA.getSize()) {
        throw std::runtime_error(
                "Dimension mismatch: rhsC must be size n.");
    }

    const Grid &grid = eta.getGrid();
    double dx = grid.dx;
    double dy = grid.dy;
    double dz = grid.dz;

    double dCoef = 0;
    switch (derivativeDirection) {
        case Axis::X: dCoef = 1.0 / (dx * dx); break;
        case Axis::Y: dCoef = 1.0 / (dy * dy); break;
        case Axis::Z: dCoef = 1.0 / (dz * dz); break;
    }
    double dt_nu_over_2 = simData.dt * simData.nu * 0.5;

    //==================================================================================================================
    // --- DEFAULT SYSTEM COEFFICIENTS (same stencil; boundaries will be then overwritten) -----------------------------
    //==================================================================================================================
    for (long i = 0; i < matA.getSize(); i++) {
        double inv_k = simData.inv_k(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                                     derivativeDirection, i);
        double beta = 1 + (dt_nu_over_2 * inv_k);
        double gamma = dt_nu_over_2 / beta;

        // --- MATRIX COEFFICIENTS -------------------------------------------------------------------------------------
        subdiag[i] = -gamma * dCoef;
        supdiag[i] = -gamma * dCoef;
        diag[i] = 1 + 2 * gamma * dCoef;

        // --- RHS COEFFICIENTS ----------------------------------------------------------------------------------------
        double eta_prec_val_m1 = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                                     derivativeDirection, i - 1);
        double eta_prec_val = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                                  derivativeDirection, i);
        double eta_prec_val_p1 = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                                     derivativeDirection, i + 1);
        double d2_eta_prec = (eta_prec_val_p1 + eta_prec_val_m1 - 2.0 * eta_prec_val) * dCoef;
        double xi_val = xi(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                           derivativeDirection, i);
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
        double physical_x = grid.to_x(iStart, eta(fieldComponent).getOffset(), eta(fieldComponent).getOffsetAxis());
        double physical_y = grid.to_y(jStart, eta(fieldComponent).getOffset(), eta(fieldComponent).getOffsetAxis());
        double physical_z = grid.to_z(kStart, eta(fieldComponent).getOffset(), eta(fieldComponent).getOffsetAxis());

        std::function<double(double, double, double, double)> bc;
        switch (fieldComponent)
        {
            case Axis::X: bc = simData.bcu; break;
            case Axis::Y: bc = simData.bcv; break;
            case Axis::Z: bc = simData.bcw; break;
        }

        switch (boundaryType)
        {
            case BoundaryType::Normal:
            {
                // Force Wall Coordinate for Normal BC (LinearSys logic)
                if (derivativeDirection == Axis::X) physical_x = 0.0;
                if (derivativeDirection == Axis::Y) physical_y = 0.0;
                if (derivativeDirection == Axis::Z) physical_z = 0.0;

                double inv_k = simData.inv_k(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                                             derivativeDirection, 0);
                double beta = 1.0 + (dt_nu_over_2 * inv_k);
                double gamma = dt_nu_over_2 / beta;

                diag.front() = 1.0 + 4.0 * gamma * dCoef;
                supdiag.front() = -4.0 / 3.0 * gamma * dCoef;

                double val_0 = xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, 0);
                double eta_0 = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, 0);
                double eta_1 = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, 1);

                rhsC.front() = val_0
                               + 4.0 / 3.0 * gamma * dCoef * eta_1
                               - 4.0 * gamma * dCoef * eta_0
                               + 8.0 / 3.0 * gamma * dCoef * (
                        bc(physical_x, physical_y, physical_z, simData.currTime - simData.dt) +
                        bc(physical_x, physical_y, physical_z, simData.currTime)
                );

                break;
            }
            case BoundaryType::Tangent:
            {
                diag.front() = 1.0;
                supdiag.front() = 0.0;

                rhsC.front() = bc(physical_x, physical_y, physical_z, simData.currTime);

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
                                      eta(fieldComponent).getOffset(),
                                      eta(fieldComponent).getOffsetAxis());
        double physical_y = grid.to_y(jStart + (derivativeDirection == Axis::Y ? endI : 0),
                                      eta(fieldComponent).getOffset(),
                                      eta(fieldComponent).getOffsetAxis());
        double physical_z = grid.to_z(kStart + (derivativeDirection == Axis::Z ? endI : 0),
                                      eta(fieldComponent).getOffset(),
                                      eta(fieldComponent).getOffsetAxis());

        std::function<double(double, double, double, double)> bc;
        switch (fieldComponent)
        {
            case Axis::X: bc = simData.bcu; break;
            case Axis::Y: bc = simData.bcv; break;
            case Axis::Z: bc = simData.bcw; break;
        }

        switch (boundaryType)
        {
            case BoundaryType::Normal:
            {
                diag.back() = 1.0;
                subdiag.back() = 0.0;

                rhsC.back() = bc(physical_x, physical_y, physical_z, simData.currTime);

                break;
            }
            case BoundaryType::Tangent:
            {
                // Enforce physical position on the actual wall
                if (derivativeDirection == Axis::X) physical_x += 0.5 * dx;
                if (derivativeDirection == Axis::Y) physical_y += 0.5 * dy;
                if (derivativeDirection == Axis::Z) physical_z += 0.5 * dz;

                double inv_k = simData.inv_k(fieldComponent).valueWithOffset(iStart, jStart, kStart,
                                                                             derivativeDirection, endI);
                double beta = 1.0 + (dt_nu_over_2 * inv_k);
                double gamma = dt_nu_over_2 / beta;

                diag.back() = 1.0 + 4.0 * gamma * dCoef;
                subdiag.back() = -4.0 / 3.0 * gamma * dCoef;

                double val_N = xi(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, endI);
                double eta_N = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, endI);
                double eta_Nm1 = eta(fieldComponent).valueWithOffset(iStart, jStart, kStart, derivativeDirection, endI - 1);

                rhsC.back() = val_N
                              + 4.0 / 3.0 * gamma * dCoef * eta_Nm1
                              - 4.0 * gamma * dCoef * eta_N
                              + 8.0 / 3.0 * gamma * dCoef * (
                                  bc(physical_x, physical_y, physical_z, simData.currTime - simData.dt) +
                                  bc(physical_x, physical_y, physical_z, simData.currTime)
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