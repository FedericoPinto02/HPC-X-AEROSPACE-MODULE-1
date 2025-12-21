#include "simulation/pressureStep.hpp"

PressureStep::PressureStep(MpiEnv &mpi, SimulationData &simData)
        : mpi(mpi), haloHandler(mpi), data_(simData) {}


void PressureStep::setup() {
    // --- Setup pressure-like fields ----------------------------------------------------------------------------------
    psi.setup(data_.grid);
    phi.setup(data_.grid);
    pcr.setup(data_.grid);
    divU.setup(data_.grid);

    // --- Assemble linear system matrices and initialize solvers for each direction -----------------------------------
    TridiagMat matrix(data_.grid->Nx);
    assembleLocalMatrix(data_.grid, Axis::X, matrix);
    solver_x = std::make_unique<SchurSolver>(mpi, Axis::X, matrix);
    solver_x->preprocess();

    matrix.resize(data_.grid->Ny);
    assembleLocalMatrix(data_.grid, Axis::Y, matrix);
    solver_y = std::make_unique<SchurSolver>(mpi, Axis::Y, matrix);
    solver_y->preprocess();

    matrix.resize(data_.grid->Nz);
    assembleLocalMatrix(data_.grid, Axis::Z, matrix);
    solver_z = std::make_unique<SchurSolver>(mpi, Axis::Z, matrix);
    solver_z->preprocess();
}


void PressureStep::run() {
    // Every pressure-like variable is solved exploting Neumann homogeneous boundary conditions

    Axis normalAxis;        // solve on the face orthogonal to normalAxis

    //==================================================================================================================
    // --- 1. Solve Psi (Direction X) ----------------------------------------------------------------------------------
    // --- when solving Psi we fill linsys with dxx derivatives
    //==================================================================================================================
    {
        haloHandler.exchange(data_.u);
        deriv.computeDivergence(data_.u, divU, data_.bcu, data_.bcv, data_.bcw, data_.currTime);

        double inv_dt = 1.0 / data_.dt;
        normalAxis = Axis::X;

        // --- Scratch linear system vectors resizing
        size_t sysDimension = data_.grid->Nx;
        rhs.resize(sysDimension);
        solution.resize(sysDimension);

        for (size_t k = 0; k < data_.grid->Nz; k++)
        {
            for (size_t j = 0; j < data_.grid->Ny; j++)
            {
                for (size_t i = 0; i < sysDimension; i++)
                {
                    rhs[i] = -divU(i, j, k) * inv_dt;
                }
                if (!data_.grid->hasMinBoundary(normalAxis))
                {
                    rhs.front() *= 0.5;
                }
                if (!data_.grid->hasMaxBoundary(normalAxis))
                {
                    rhs.back() *= 0.5;
                }

                solver_x->solve(rhs, solution);

                for (size_t i = 0; i < sysDimension; i++)
                {
                    psi(i, j, k) = solution[i];
                }
            }
        }
    }


    //==================================================================================================================
    // --- 2. Solve Phi (Direction Y) ----------------------------------------------------------------------------------
    // --- when solving Phi we fill linsys with dyy derivatives
    //==================================================================================================================
    {
        normalAxis = Axis::Y;

        // --- Scratch linear system vectors resizing
        size_t sysDimension = data_.grid->Ny;
        rhs.resize(sysDimension);
        solution.resize(sysDimension);

        for (size_t k = 0; k < data_.grid->Nz; k++)
        {
            for (size_t i = 0; i < data_.grid->Nx; i++)
            {
                for (size_t j = 0; j < sysDimension; j++)
                {
                    rhs[j] = psi(i, j, k);
                }
                if (!data_.grid->hasMinBoundary(normalAxis))
                {
                    rhs.front() *= 0.5;
                }
                if (!data_.grid->hasMaxBoundary(normalAxis))
                {
                    rhs.back() *= 0.5;
                }

                solver_y->solve(rhs, solution);

                for (size_t j = 0; j < sysDimension; j++)
                {
                    phi(i, j, k) = solution[j];
                }
            }
        }
    }


    //==================================================================================================================
    // --- 3. Solve PCR (second phi) (Direction Z) ---------------------------------------------------------------------
    // --- when solving Pcr we fill linsys with dzz derivatives
    //==================================================================================================================
    {
        normalAxis = Axis::Z;   // x = pcr, rhs = phi

        // --- Scratch linear system vectors resizing
        size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve
        rhs.resize(sysDimension);
        solution.resize(sysDimension);

        for (size_t j = 0; j < data_.grid->Ny; j++)
        {
            for (size_t i = 0; i < data_.grid->Nx; i++)
            {
                for (size_t k = 0; k < sysDimension; k++)
                {
                    rhs[k] = phi(i, j, k);
                }
                if (!data_.grid->hasMinBoundary(normalAxis))
                {
                    rhs.front() *= 0.5;
                }
                if (!data_.grid->hasMaxBoundary(normalAxis))
                {
                    rhs.back() *= 0.5;
                }

                solver_z->solve(rhs, solution);

                for (size_t k = 0; k < sysDimension; k++)
                {
                    pcr(i, j, k) = solution[k];
                }
            }
        }
    }


    //==================================================================================================================
    // --- 4. Update pressure with rotational correction, derive pressure predictor ------------------------------------
    //==================================================================================================================
    double chi = 0.0;
    double Re = 1;
    double nu = data_.nu;
    double factor = chi * nu;

    for (size_t k = 0; k < data_.grid->Nz; k++)
    {
        for (size_t j = 0; j < data_.grid->Ny; j++)
        {
            for (size_t i = 0; i < data_.grid->Nx; i++)
            {
                data_.p(i, j, k) += pcr(i, j, k) - factor * divU(i, j, k);
                data_.predictor(i, j, k) = data_.p(i, j, k) + pcr(i, j, k);
            }
        }
    }
}


void PressureStep::assembleLocalMatrix(
        const GridPtr &grid,
        const Axis direction,
        TridiagMat &matrix
) {
    std::vector<double> &diag = matrix.getDiag(0);
    std::vector<double> &subdiag = matrix.getDiag(-1);
    std::vector<double> &supdiag = matrix.getDiag(1);

    double delta = 0;
    switch (direction) {
        case Axis::X: delta = grid->dx; break;
        case Axis::Y: delta = grid->dy; break;
        case Axis::Z: delta = grid->dz; break;
    }
    double inv_delta = 1.0 / delta;
    double inv_delta2 = inv_delta * inv_delta;

    std::fill(subdiag.begin(), subdiag.end(), -inv_delta2);
    std::fill(supdiag.begin(), supdiag.end(), -inv_delta2);
    std::fill(diag.begin(), diag.end(), 1.0 + 2.0 * inv_delta2);

    // Boundary conditions (Neumann homogeneous)
    /*
     * (a0) b0 c0
     *      a1 b1 c1
     *         a2 b2 c2
     *            a3 b3 (c3)
     */
    // --- LEFT INTERFACE ---
    if (grid->hasMinBoundary(direction)) {
        // physical
        supdiag.front() = -2.0 * inv_delta2;
    }
    else {
        // internal (half contribution)
        diag.front() *= 0.5;
    }
    // --- RIGHT INTERFACE ---
    if (grid->hasMaxBoundary(direction)) {
        // physical
        diag.back() = 1.0 + inv_delta2;
    }
    else {
        // internal (half contribution)
        diag.back() *= 0.5;
    }
}
