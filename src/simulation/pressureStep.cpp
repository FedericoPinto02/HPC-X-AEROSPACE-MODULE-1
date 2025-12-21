#include "simulation/pressureStep.hpp"

PressureStep::PressureStep(MpiEnv &mpi, SimulationData &simData)
        : mpi(mpi), haloHandler(mpi), data_(simData) {
    psi.setup(data_.grid);
    phi.setup(data_.grid);
    pcr.setup(data_.grid);
    divU.setup(data_.grid);
}


void PressureStep::run() {
    // Linear system variables definition: progressive resizing instead of per-sweep reallocation
    TridiagMat matrix;
    std::vector<double> rhs, solution;


    Axis normalAxis;        // solve on the face orthogonal to normalAxis
    // Every pressure-like variable is solved exploting Neumann homogeneous boundary conditions


    // ------------------------------------------
    // Solve Psi --------------------------------
    // ------------------------------------------
    {
        double inv_dt = 1.0 / data_.dt;

        normalAxis = Axis::X;   // x = psi, rhs = divU/dt
        size_t sysDimension = data_.grid->Nx;

        haloHandler.exchange(data_.u);
        deriv.computeDivergence(data_.u, divU, data_.bcu, data_.bcv, data_.bcw, data_.currTime);

        // --- Linear system and solver definition
        matrix.resize(sysDimension);
        rhs.resize(sysDimension);
        solution.resize(sysDimension);

        // when solving Psi we fill linsys with dxx derivatives
        assembleLocalMatrix(data_.grid, normalAxis, matrix);
        SchurSolver solver(mpi, normalAxis, matrix);
        solver.preprocess();

        for (size_t k = 0; k < data_.grid->Nz; k++)
        {
            for (size_t j = 0; j < data_.grid->Ny; j++)
            {
                for (size_t i = 0; i < sysDimension; i++)
                {
                    rhs[i] = -divU(i, j, k) * inv_dt;
                }
                if (!data_.grid->hasMinBoundary(Axis::X))
                {
                    rhs.front() *= 0.5;
                }
                if (!data_.grid->hasMaxBoundary(Axis::X))
                {
                    rhs.back() *= 0.5;
                }

                solver.solve(rhs, solution);

                for (size_t i = 0; i < sysDimension; i++)
                {
                    psi(i, j, k) = solution[i];
                }
            }
        }
    }


    // ------------------------------------------
    // Solve Phi --------------------------------
    // ------------------------------------------
    {
        normalAxis = Axis::Y;   // x = phi, rhs = psi
        size_t sysDimension = data_.grid->Ny;

        // --- Linear system and solver definition
        matrix.resize(sysDimension);
        rhs.resize(sysDimension);
        solution.resize(sysDimension);

        // when solving Phi we fill linsys with dyy derivatives
        assembleLocalMatrix(data_.grid, normalAxis, matrix);
        SchurSolver solver(mpi, normalAxis, matrix);
        solver.preprocess();

        for (size_t k = 0; k < data_.grid->Nz; k++)
        {
            for (size_t i = 0; i < data_.grid->Nx; i++)
            {
                for (size_t j = 0; j < sysDimension; j++)
                {
                    rhs[j] = psi(i, j, k);
                }
                if (!data_.grid->hasMinBoundary(Axis::Y))
                {
                    rhs.front() *= 0.5;
                }
                if (!data_.grid->hasMaxBoundary(Axis::Y))
                {
                    rhs.back() *= 0.5;
                }

                solver.solve(rhs, solution);

                for (size_t j = 0; j < sysDimension; j++)
                {
                    phi(i, j, k) = solution[j];
                }
            }
        }
    }

    // ------------------------------------------
    // Solve Pcr --------------------------------
    // ------------------------------------------
    {
        normalAxis = Axis::Z;   // x = pcr, rhs = phi
        size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve

        // --- Linear system and solver definition
        matrix.resize(sysDimension);
        rhs.resize(sysDimension);
        solution.resize(sysDimension);

        // when solving Pcr we fill linsys with dzz derivatives
        assembleLocalMatrix(data_.grid, normalAxis, matrix);
        SchurSolver solver(mpi, normalAxis, matrix);
        solver.preprocess();

        for (size_t j = 0; j < data_.grid->Ny; j++)
        {
            for (size_t i = 0; i < data_.grid->Nx; i++)
            {
                for (size_t k = 0; k < sysDimension; k++)
                {
                    rhs[k] = phi(i, j, k);
                }
                if (!data_.grid->hasMinBoundary(Axis::Z))
                {
                    rhs.front() *= 0.5;
                }
                if (!data_.grid->hasMaxBoundary(Axis::Z))
                {
                    rhs.back() *= 0.5;
                }

                solver.solve(rhs, solution);

                for (size_t k = 0; k < sysDimension; k++)
                {
                    pcr(i, j, k) = solution[k];
                }
            }
        }
    }

    // Add pressure corrector contribution and rotational correction
    // Write pressure predictor
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

    // Boundary conditions
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
