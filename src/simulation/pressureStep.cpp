#include <algorithm>
#include <cmath>
#include <vector>
#include <core/Fields.hpp>
#include <numerics/derivatives.hpp>
#include <numerics/LinearSys.hpp>
#include <simulation/pressureStep.hpp>



PressureStep::PressureStep(SimulationData& simData) : data_(simData)
{
    initializeWorkspaceFields(data_.gridPtr);
}



void PressureStep::initializeWorkspaceFields(std::shared_ptr<const Grid> gridPtr)
{
    std::vector<Field::Scalar> zeros(gridPtr->size(), 0.0);
    psi.setup(gridPtr, zeros);
    phi.setup(gridPtr, zeros);
    pcr.setup(gridPtr, zeros);
    divU.setup(gridPtr, zeros);
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
    size_t nSystem = data_.gridPtr->Ny * data_.gridPtr->Nz; // number of linear systems to solve
    size_t sysDimension = data_.gridPtr->Nx; // dimension of linear system to solve
    // when solving Psi we fill linsys with dxx derivatives

    LinearSys mySystem(sysDimension, BoundaryType::Normal);
    std::vector<double> rhs(sysDimension);
    rhs.front() = 0;
    rhs.back() = 0;

    
    Derivatives deriv;
    deriv.computeDivergence(data_.u, divU);
    double inv_dt = 1.0 / data_.dt;

    // iStart = 0;   
    for (size_t j = 0; j < data_.gridPtr->Ny; j++)
    {
        for (size_t k = 0; k < data_.gridPtr->Nz; k++)
        {
            // jStart = j;
            // kStart = k;

            for (size_t i = 1; i < sysDimension-1; i++)
            {
                rhs[i] = divU(i,j,k) * inv_dt;
            }

            mySystem.setRhs(rhs);

            mySystem.fillSystemPressure(data_.p, normalAxis);

            mySystem.ThomaSolver();

            std::vector<Field::Scalar> solution = mySystem.getSolution();
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
    size_t nSystem = data_.gridPtr->Nx * data_.gridPtr->Nz; // number of linear systems to solve
    size_t sysDimension = data_.gridPtr->Ny; // dimension of linear system to solve
    // when solving Phi we fill linsys with dyy derivatives

    LinearSys mySystem(sysDimension, BoundaryType::Normal);
    std::vector<double> rhs(sysDimension);      

    // jStart = 0;
    for (size_t i = 0; i < data_.gridPtr->Nx; i++)
    {
        for (size_t k = 0; k < data_.gridPtr->Nz; k++)
        {
            // iStart = i;
            // kStart = k;

            for (size_t j = 1; j < sysDimension-1; j++)
            {
                rhs[j] = psi(i,j,k);
            }

            // Double check on boundary conditions
            rhs.front() = 0.0;
            rhs.back() = 0.0;

            mySystem.setRhs(rhs);

            mySystem.fillSystemPressure(data_.p, normalAxis);

            mySystem.ThomaSolver();

            std::vector<Field::Scalar> solution = mySystem.getSolution();
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
    size_t nSystem = data_.gridPtr->Nx * data_.gridPtr->Ny; // number of linear systems to solve
    size_t sysDimension = data_.gridPtr->Nz; // dimension of linear system to solve
    // when solving Pcr we fill linsys with dzz derivatives

    LinearSys mySystem(sysDimension, BoundaryType::Normal);
    std::vector<double> rhs(sysDimension);      

    // kStart = 0;
    for (size_t j = 0; j < data_.gridPtr->Ny; j++)
    {
        for (size_t i = 0; i < data_.gridPtr->Nx; i++)
        {
            // jStart = j;
            // iStart = i;

            for (size_t k = 1; k < sysDimension-1; k++)
            {
                rhs[k] = phi(i,j,k);
            }

            // Double check on boundary conditions
            rhs.front() = 0.0;
            rhs.back() = 0.0;

            mySystem.setRhs(rhs);

            mySystem.fillSystemPressure(data_.p, normalAxis);

            mySystem.ThomaSolver();

            std::vector<Field::Scalar> solution = mySystem.getSolution();
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


