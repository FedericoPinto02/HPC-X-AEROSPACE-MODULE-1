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



PressureStep::PressureStep(SimulationData& simData, ParallelizationSettings& parallel) : data_(simData), parallel_(parallel)
{
    initializeWorkspaceFields();
}



void PressureStep::initializeWorkspaceFields()
{
    psi.setup(data_.grid);
    phi.setup(data_.grid);
    pcr.setup(data_.grid);
    divU.setup(data_.grid);
}


void PressureStep::run()
{   
   
    // size_t iStart, jStart, kStart;   // these will be useful with parallelization
    // solve on the face orthogonal to normalAxis
    Axis normalAxis;
    // Every pressure-like variable is solved exploting Neumann homogeneous boundary conditions
    
    
    // ------------------------------------------
    // Solve Psi --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::X;   // x = psi, rhs = divU/dt
    size_t nSystem = data_.grid->Ny * data_.grid->Nz; // number of linear systems to solve
    size_t sysDimension = data_.grid->Nx; // dimension of linear system to solve
    // when solving Psi we fill linsys with dxx derivatives

    LinearSys mySystem(sysDimension, BoundaryType::Normal);
    std::vector<double> rhs(sysDimension);

    
    Derivatives deriv;
    deriv.computeDivergence(data_.u, divU, data_);
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

            mySystem.setRhs(rhs);

            mySystem.fillSystemPressure(data_.p, normalAxis);

            // mySystem.ThomaSolver();

            // std::vector<Field::Scalar> solution = mySystem.getSolution();
            std::vector<Field::Scalar> solution = solveSystem(mySystem, BoundaryType::Normal);
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
    size_t nSystem = data_.grid->Nx * data_.grid->Nz; // number of linear systems to solve
    size_t sysDimension = data_.grid->Ny; // dimension of linear system to solve
    // when solving Phi we fill linsys with dyy derivatives

    LinearSys mySystem(sysDimension, BoundaryType::Normal);
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

            mySystem.setRhs(rhs);

            mySystem.fillSystemPressure(data_.p, normalAxis);

            // mySystem.ThomaSolver();
            
            // std::vector<Field::Scalar> solution = mySystem.getSolution();
            std::vector<Field::Scalar> solution = solveSystem(mySystem, BoundaryType::Normal);
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
    size_t nSystem = data_.grid->Nx * data_.grid->Ny; // number of linear systems to solve
    size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve
    // when solving Pcr we fill linsys with dzz derivatives

    LinearSys mySystem(sysDimension, BoundaryType::Normal);
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

            mySystem.setRhs(rhs);

            mySystem.fillSystemPressure(data_.p, normalAxis);

            // mySystem.ThomaSolver();

            // std::vector<Field::Scalar> solution = mySystem.getSolution();
            std::vector<Field::Scalar> solution = solveSystem(mySystem, BoundaryType::Normal);
            for (size_t k = 0; k < sysDimension; k++)
            {
                pcr(i, j, k) = 1.0 * solution[k];
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

    for (size_t k = 0; k < data_.grid->Nz; k++) {
    for (size_t j = 0; j < data_.grid->Ny; j++) {
        for (size_t i = 0; i < data_.grid->Nx; i++) {
            
            data_.p(i, j, k) += pcr(i, j, k) - factor * divU(i, j, k);
            data_.predictor(i,j,k) = data_.p(i, j, k) + pcr(i, j, k);
        }
    }
    
   
}
}


