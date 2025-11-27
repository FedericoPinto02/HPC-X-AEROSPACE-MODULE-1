#include <algorithm>
#include <cmath>
#include <vector>
#include <core/Fields.hpp>
#include <numerics/derivatives.hpp>
#include <numerics/LinearSys.hpp>
#include "numerics/SchurSequentialSolver.hpp"
#include <simulation/pressureStep.hpp>


std::vector<double> PressureStep::solveSystem(LinearSys& sys, BoundaryType bType) {

    // 1. If we want to use the classic solver (debug or actual P=1)
    if (parallel_.schurDomains <= 1) {
        sys.ThomaSolver();
        return sys.getSolution();
    }
    else 
    {
        // 2. Sequential Schur Logic
        SchurSequentialSolver schur(sys.getMatrix().getSize(), parallel_.schurDomains, bType);
        schur.PreProcess(sys.getMatrix());
        return schur.solve(sys.getRhs());
    }
}


PressureStep::PressureStep(SimulationData &simData, ParallelizationSettings &parallel)
        : data_(simData), parallel_(parallel),
          linearSys_xLine(simData.grid->Nx, BoundaryType::Normal),
          linearSys_yLine(simData.grid->Ny, BoundaryType::Normal),
          linearSys_zLine(simData.grid->Nz, BoundaryType::Normal)
{
    psi.setup(data_.grid);
    phi.setup(data_.grid);
    pcr.setup(data_.grid);
    divU.setup(data_.grid);
    assembleLineMatrix(data_.grid->dx, linearSys_xLine);
    assembleLineMatrix(data_.grid->dy, linearSys_yLine);
    assembleLineMatrix(data_.grid->dz, linearSys_zLine);
}

void PressureStep::assembleLineMatrix(double delta, LinearSys& sys)
{
    double inv_delta = 1.0 / (delta * delta);
    TridiagMat &mat = sys.getMatrix();
    auto &diag = mat.getDiag(0), &subdiag = mat.getDiag(-1), &supdiag = mat.getDiag(1);
    std::fill(subdiag.begin(), subdiag.end(), -inv_delta);
    std::fill(supdiag.begin(), supdiag.end(), -inv_delta);
    std::fill(diag.begin(), diag.end(), 1.0 + 2.0 * inv_delta);

    // Boundary conditions
    supdiag.front() = -2.0 * inv_delta;
    diag.back() = 1.0 + inv_delta;
}


void PressureStep::run()
{   
    // size_t iStart, jStart, kStart;   // these will be useful with parallelization
    // Every pressure-like variable is solved exploiting Neumann homogeneous boundary conditions
    
    
    // ------------------------------------------
    // Solve Psi --------------------------------
    // ------------------------------------------
    {
    size_t sysDimension = data_.grid->Nx; // dimension of linear system to solve
    // when solving Psi we fill linsys with dxx derivatives

    std::vector<double> rhs(sysDimension);

    
    Derivatives deriv;
    deriv.computeDivergence(data_.u, divU);
    double inv_dt = 1.0 / data_.dt;

    // iStart = 0;   
    for (size_t j = 0; j < data_.grid->Ny; j++)
    {
        for (size_t k = 0; k < data_.grid->Nz; k++)
        {
            // jStart = j;
            // kStart = k;

            for (size_t i = 0; i < sysDimension; i++)
            {
                rhs[i] =  - divU(i,j,k) * inv_dt;
            }

            linearSys_xLine.setRhs(rhs);

            std::vector<Field::Scalar> solution = solveSystem(linearSys_xLine, BoundaryType::Normal);
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
    size_t sysDimension = data_.grid->Ny; // dimension of linear system to solve
    // when solving Phi we fill linsys with dyy derivatives

    std::vector<double> rhs(sysDimension);

    // jStart = 0;
    for (size_t i = 0; i < data_.grid->Nx; i++)
    {
        for (size_t k = 0; k < data_.grid->Nz; k++)
        {
            // iStart = i;
            // kStart = k;

            for (size_t j = 0; j < sysDimension; j++)
            {
                rhs[j] = psi(i,j,k);
            }

            linearSys_yLine.setRhs(rhs);

            std::vector<Field::Scalar> solution = solveSystem(linearSys_yLine, BoundaryType::Normal);
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
    size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve
    // when solving Pcr we fill linsys with dzz derivatives

    std::vector<double> rhs(sysDimension);

    // kStart = 0;
    for (size_t j = 0; j < data_.grid->Ny; j++)
    {
        for (size_t i = 0; i < data_.grid->Nx; i++)
        {
            // jStart = j;
            // iStart = i;

            for (size_t k = 0; k < sysDimension; k++)
            {
                rhs[k] = phi(i,j,k);
            }

            linearSys_zLine.setRhs(rhs);

            std::vector<Field::Scalar> solution = solveSystem(linearSys_zLine, BoundaryType::Normal);
            for (size_t k = 0; k < sysDimension; k++)
            {
                pcr(i, j, k) = solution[k];
            }
           
        }
    }
    }

    // Add pressure corrector contribution
    data_.p.add( pcr );

   
}


