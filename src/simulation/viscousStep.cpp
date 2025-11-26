#include <algorithm>
#include <cmath>
#include <core/Fields.hpp>
#include <vector>
#include <numerics/derivatives.hpp>
#include <numerics/LinearSys.hpp>
#include "numerics/SchurSequentialSolver.hpp"
#include <simulation/viscousStep.hpp>


std::vector<double> ViscousStep::solveSystem(LinearSys& sys, BoundaryType bType) {

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


ViscousStep::ViscousStep(SimulationData& simData, ParallelizationSettings& parallel) : data_(simData), parallel_(parallel)
{
    initializeWorkspaceFields();
}



void ViscousStep::run()
{
    computeG();
    computeXi();
    closeViscousStep();
}

void ViscousStep::initializeWorkspaceFields()
{
    g.setup(data_.grid);
    gradP.setup(data_.grid);
    dxxEta.setup(data_.grid);
    dyyZeta.setup(data_.grid);
    dzzU.setup(data_.grid);
    xi.setup(data_.grid);
}

void ViscousStep::computeG()
{
    // Ingredients list
    auto& eta = data_.eta;
    auto& zeta = data_.zeta;
    auto& u = data_.u;
    auto& predictor = data_.predictor;
    auto& k_data = data_.k.getData(); // k cant be 0!!!
    double nu_val = data_.nu;
    
    Derivatives derive;
    derive.computeGradient(predictor, gradP);
    

    // Recepie
    // g = f  -grad(p)  -nu/k *u * 0.5  +nu*(dxx eta + dyy zeta + dzz u) * 0.5

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {

        for ( size_t k=0; k < data_.grid->Nz; k++)
            for ( size_t j=0; j < data_.grid->Ny; j++)
                for ( size_t i=0; i < data_.grid->Nx; i++)
                {
                    dxxEta(axis, i, j, k) = derive.Dxx_local(eta(axis), i, j, k);
                    dyyZeta(axis, i, j, k) = derive.Dyy_local(zeta(axis), i, j, k);
                    dzzU(axis, i, j, k) = derive.Dzz_local(u(axis), i, j, k);
                }

        auto& f_data = data_.f(axis).getData();
        auto& u_data = data_.u(axis).getData();
        auto& gradP_data = gradP(axis).getData();
        auto& dxx_data = dxxEta(axis).getData();
        auto& dyy_data = dyyZeta(axis).getData();
        auto& dzz_data = dzzU(axis).getData();
        auto& g_data = g(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++)
        {
            g_data[i] = f_data[i] - gradP_data[i] - nu_val * u_data[i] / k_data[i] * 0.5 + nu_val * (dxx_data[i] + dyy_data[i] + dzz_data[i]) * 0.5;
        }

        const auto &grid = *data_.grid;
        std::vector<double> mx(grid.Nx), my(grid.Ny), mz(grid.Nz);

        for (Axis axis : {Axis::X, Axis::Y, Axis::Z})
        {
            for (size_t i = 0; i < grid.Nx; i++)
                mx[i] = (i == 0 || (axis == Axis::X && i == grid.Nx - 1)) ? 0.0 : 1.0;

            for (size_t j = 0; j < grid.Ny; j++)
                my[j] = (j == 0 || (axis == Axis::Y && j == grid.Ny - 1)) ? 0.0 : 1.0;

            for (size_t k = 0; k < grid.Nz; k++)
                mz[k] = (k == 0 || (axis == Axis::Z && k == grid.Nz - 1)) ? 0.0 : 1.0;

                
            for (size_t k = 0; k < grid.Nz; k++)
            {
                const double mk = mz[k];

                for (size_t j = 0; j < grid.Ny; j++)
                {
                    const double mj = mk * my[j];

                    for (size_t i = 0; i < grid.Nx; i++)
                    {
                        g(axis, i, j, k) *= mj * mx[i];
                    }
                }
            }
        }
    }
}

void ViscousStep::computeXi()
{
    // Ingredients list
    auto& u = data_.u;
    auto& k_data = data_.k.getData();  // k cant be 0!!!
    double nu_val = data_.nu;
    double dt_val = data_.dt;
   
    // Recepie
    // beta = 1+ dt*nu /2/k
    // xi = u + dt/beta * g

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {

        auto& u_data = data_.u(axis).getData();
        auto& g_data = g(axis).getData();
        auto& xi_data = xi(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++)
        {
            xi_data[i] = u_data[i]
                        + dt_val * g_data[i]
                        / (1 + dt_val * nu_val * 0.5 / k_data[i]);
        }
    }
}

void ViscousStep::closeViscousStep()
{   
    const double dx = data_.grid->dx;
    const double dy = data_.grid->dy;
    const double dz = data_.grid->dz;
    const double t = data_.currTime;
    double porosity, beta, gamma, mul, deriv_u, deriv_v, deriv_w;
    size_t iStart, jStart, kStart;
    // solve on the face orthogonal to normalAxis
    Axis normalAxis;
    

    // ------------------------------------------
    // Solve Eta --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::X;
    size_t nSystem = data_.grid->Ny * data_.grid->Nz; // number of linear systems to solve
    size_t sysDimension = data_.grid->Nx; // dimension of linear system to solve
    size_t j, k;
    // when solving Eta we fill linsys with dxx derivatives
    // Eta.u is then solved exploiting normal Dirichlet boundary conditions
    // Eta.v and Eta.w are solved exploiting tangent Dirichlet boundary conditions
    // Differences between normal and tangent are given by staggered grid.

    LinearSys mySystem_u(sysDimension, BoundaryType::Normal);
    LinearSys mySystem_v(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_w(sysDimension, BoundaryType::Tangent);
    std::vector<double> rhs_u(sysDimension);
    std::vector<double> rhs_v(sysDimension);
    std::vector<double> rhs_w(sysDimension);

    std::vector<double> unknown_u;
    unknown_u.reserve(sysDimension);
    std::vector<double> unknown_v;
    unknown_v.reserve(sysDimension);
    std::vector<double> unknown_w;
    unknown_w.reserve(sysDimension);

    rhs_u.front() = 0;
    rhs_u.back() = 0;
    rhs_v.front() = 0;
    rhs_v.back() = 0;
    rhs_w.front() = 0;
    rhs_w.back() = 0;

    iStart = 0;
    mul = 1.0 / (data_.grid->dx * data_.grid->dx);

    for (j = 1; j < data_.grid->Ny-1; j++)
    {
        for (k = 1; k < data_.grid->Nz-1; k++)
        {
            jStart = j;
            kStart = k;

            for (size_t i = 1; i < sysDimension-1; i++)
            {
                porosity = data_.k(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                deriv_u = (data_.eta(Axis::X, i + 1, j, k) + data_.eta(Axis::X, i - 1, j, k) - 2.0 * data_.eta(Axis::X, i, j, k))*mul;
                deriv_v = (data_.eta(Axis::Y, i + 1, j, k) + data_.eta(Axis::Y, i - 1, j, k) - 2.0 * data_.eta(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.eta(Axis::Z, i + 1, j, k) + data_.eta(Axis::Z, i - 1, j, k) - 2.0 * data_.eta(Axis::Z, i, j, k))*mul;

                rhs_u[i] = xi(Axis::X, i,j,k) - gamma * deriv_u;
                rhs_v[i] = xi(Axis::Y, i,j,k) - gamma * deriv_v;
                rhs_w[i] = xi(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);
            
            mySystem_u.fillSystemVelocity(data_, data_.eta, xi, Axis::X, Axis::X, iStart, jStart, kStart);
            // mySystem_u.ThomaSolver();
            // std::vector<double> unknown_u = mySystem_u.getSolution();
            unknown_u = solveSystem(mySystem_u, BoundaryType::Normal);


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.eta, xi, Axis::Y, Axis::X, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.eta, xi, Axis::Z, Axis::X, iStart, jStart, kStart);
            // mySystem_w.ThomaSolver();
            // std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);
            
            for (size_t i = 0; i < sysDimension; i++)
            {
                data_.eta(Axis::X, i, j, k) = unknown_u[i];
                data_.eta(Axis::Y, i, j, k) = unknown_v[i];
                data_.eta(Axis::Z, i, j, k) = unknown_w[i];
            }

        }
    }
    // Special case, k = 1/2 (=0)
    // Here Eta.w is computed with tangent b.c while Eta.u and Eta.v come directly from boundary conditions.
    k = 0;
    iStart = 0;
    for (j = 1; j < data_.grid->Ny - 1; j++)
    {
        jStart = j;
        kStart = k;

        for (size_t i = 1; i < sysDimension - 1; i++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_w = (data_.eta(Axis::Z, i + 1, j, k) + data_.eta(Axis::Z, i - 1, j, k) - 2.0 * data_.eta(Axis::Z, i, j, k)) * mul;
            rhs_w[i] = xi(Axis::Z, i, j, k) - gamma * deriv_w;
        }

        mySystem_w.setRhs(rhs_w);

        mySystem_w.fillSystemVelocity(data_, data_.eta, xi, Axis::Z, Axis::X, iStart, jStart, kStart);
        unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);
        // mySystem_w.ThomaSolver();
        // std::vector<double> unknown_w = mySystem_w.getSolution();

        for (size_t i = 0; i < sysDimension; i++)
        {
            data_.eta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
            data_.eta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
            data_.eta(Axis::Z, i, j, k) = unknown_w[i];
        }
    }
    // Special case, k = Nz-1
    // Here Eta.u, Eta.v are computed with tangent b.c while , Eta.w comes directly from boundary conditions.
    k = data_.grid->Nz - 1;
    for (j = 1; j < data_.grid->Ny - 1; j++)
    {
        for (size_t i = 1; i < sysDimension - 1; i++)
        {
            jStart = j;
            kStart = k;
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_u = (data_.eta(Axis::X, i + 1, j, k) + data_.eta(Axis::X, i - 1, j, k) - 2.0 * data_.eta(Axis::X, i, j, k))*mul;
            rhs_u[i] = xi(Axis::X, i,j,k) - gamma * deriv_u;

            deriv_v = (data_.eta(Axis::Y, i + 1, j, k) + data_.eta(Axis::Y, i - 1, j, k) - 2.0 * data_.eta(Axis::Y, i, j, k))*mul;
            rhs_v[i] = xi(Axis::Y, i,j,k) - gamma * deriv_v;
        }
            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_, data_.eta, xi, Axis::X, Axis::X, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_u = solveSystem(mySystem_u, BoundaryType::Normal);

            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.eta, xi, Axis::Y, Axis::X, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);

        for (size_t i = 0; i < sysDimension; i++)
        {
            data_.eta(Axis::X, i, j, k) = unknown_u[i];
            data_.eta(Axis::Y, i, j, k) = unknown_v[i];
            data_.eta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
        }
    }

    // Special case, j = 1/2 (=0)
    // Here Eta.v is computed with tangent b.c while Eta.u and Eta.w come directly from boundary conditions.
    j = 0;
    for (k = 1; k < data_.grid->Nz - 1; k++)
    {
        for (size_t i = 1; i < sysDimension - 1; i++)
        {
            jStart = j;
            kStart = k;
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_v = (data_.eta(Axis::Y, i + 1, j, k) + data_.eta(Axis::Y, i - 1, j, k) - 2.0 * data_.eta(Axis::Y, i, j, k)) * mul;
            rhs_v[i] = xi(Axis::Y, i, j, k) - gamma * deriv_v;
        }

        mySystem_v.setRhs(rhs_v);

        mySystem_v.fillSystemVelocity(data_, xi, data_.eta, Axis::Y, Axis::X, iStart, jStart, kStart);
        //mySystem_v.ThomaSolver();
        //std::vector<double> unknown_v = mySystem_v.getSolution();
        unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);

        for (size_t i = 0; i < sysDimension; i++)
        {
            data_.eta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
            data_.eta(Axis::Y, i, j, k) = unknown_v[i];
            data_.eta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
        }
    }
    // Special case, j = Ny-1
    // Here Eta.u, Eta.w are computed with tangent b.c while Eta.v comes directly from boundary conditions.
    j = data_.grid->Ny - 1;
    for (k = 1; k < data_.grid->Nz - 1; k++)
    {
        for (size_t i = 1; i < sysDimension - 1; i++)
        {
            jStart = j;
            kStart = k;
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_u = (data_.eta(Axis::X, i + 1, j, k) + data_.eta(Axis::X, i - 1, j, k) - 2.0 * data_.eta(Axis::X, i, j, k))*mul;
            rhs_u[i] = xi(Axis::X, i,j,k) - gamma * deriv_u;

            deriv_w = (data_.eta(Axis::Z, i + 1, j, k) + data_.eta(Axis::Z, i - 1, j, k) - 2.0 * data_.eta(Axis::Z, i, j, k))*mul;
            rhs_w[i] = xi(Axis::Z, i,j,k) - gamma * deriv_w;
        }
            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_, data_.eta, xi, Axis::X, Axis::X, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_u = solveSystem(mySystem_u, BoundaryType::Normal);

            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.eta, xi, Axis::Z, Axis::X, iStart, jStart, kStart);
            //mySystem_w.ThomaSolver();
            //std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);

        for (size_t i = 0; i < sysDimension; i++)
        {
            data_.eta(Axis::X, i, j, k) = unknown_u[i];
            data_.eta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
            data_.eta(Axis::Z, i, j, k) = unknown_w[i];
        }
    }
    // --- Handle Corner Lines (where j and k are boundaries) ---
    // For these lines, no solver is run. We apply boundary conditions directly.
    // Note: sysDimension == data_.grid->Nx here.

    // Special case: Corner j = 0, k = 0
    j = 0;
    k = 0; 
    for (size_t i = 0; i < sysDimension; i++)
    {
        data_.eta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
        data_.eta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
        data_.eta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
    }

    // Special case: Corner j = Ny-1, k = Nz-1
    j = data_.grid->Ny - 1;
    k = data_.grid->Nz - 1;
    for (size_t i = 1; i < sysDimension - 1; i++)
    {
        jStart = j;
        kStart = k;
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_u = (data_.eta(Axis::X, i + 1, j, k) + data_.eta(Axis::X, i - 1, j, k) - 2.0 * data_.eta(Axis::X, i, j, k)) * mul;
        rhs_u[i] = xi(Axis::X, i, j, k) - gamma * deriv_u;
    }
    mySystem_u.setRhs(rhs_u);
    mySystem_u.fillSystemVelocity(data_, data_.eta, xi, Axis::X, Axis::X, iStart, jStart, kStart);
    unknown_u = solveSystem(mySystem_u, BoundaryType::Normal);

    for (size_t i = 0; i < sysDimension; i++)
    {
        data_.eta(Axis::X, i, j, k) = unknown_u[i];
        data_.eta(Axis::Y, i, j, k) = data_.bcv(i * dx, (j + 0.5) * dy, k * dz, t);
        data_.eta(Axis::Z, i, j, k) = data_.bcw(i * dx, j * dy, (k + 0.5) * dz, t);
    }

    // Special case: Corner j = 0, k = Nz-1
    j = 0;
    k = data_.grid->Nz - 1;
    for (size_t i = 1; i < sysDimension - 1; i++)
    {
        jStart = j;
        kStart = k;
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_v = (data_.eta(Axis::Y, i + 1, j, k) + data_.eta(Axis::Y, i - 1, j, k) - 2.0 * data_.eta(Axis::Y, i, j, k)) * mul;
        rhs_v[i] = xi(Axis::Y, i, j, k) - gamma * deriv_v;
    }
    mySystem_v.setRhs(rhs_v);
    mySystem_v.fillSystemVelocity(data_, data_.eta, xi, Axis::Y, Axis::X, iStart, jStart, kStart);
    unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);

    for (size_t i = 0; i < sysDimension; i++)
    {
        data_.eta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
        data_.eta(Axis::Y, i, j, k) = unknown_v[i];
        data_.eta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
    }

    // Special case: Corner j = Ny-1, k = 0
    j = data_.grid->Ny - 1;
    k = 0;
    for (size_t i = 1; i < sysDimension - 1; i++)
    {
        jStart = j;
        kStart = k;
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_w = (data_.eta(Axis::Z, i + 1, j, k) + data_.eta(Axis::Z, i - 1, j, k) - 2.0 * data_.eta(Axis::Z, i, j, k)) * mul;
        rhs_w[i] = xi(Axis::Z, i, j, k) - gamma * deriv_w;
    }
    mySystem_w.setRhs(rhs_w);
    mySystem_w.fillSystemVelocity(data_, data_.eta, xi, Axis::Z, Axis::X, iStart, jStart, kStart);
    unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);

    for (size_t i = 0; i < sysDimension; i++)
    {
        data_.eta(Axis::X, i, j, k) = data_.bcu((i + 0.5) * dx, j * dy, k * dz, t);
        data_.eta(Axis::Y, i, j, k) = data_.bcv(i * dx, (j + 0.5) * dy, k * dz, t);
        data_.eta(Axis::Z, i, j, k) = unknown_w[i];
    }
    }
   



    // ------------------------------------------
    // Solve Zeta --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::Y;
    size_t nSystem = data_.grid->Nx * data_.grid->Nz; // number of linear systems to solve
    size_t sysDimension = data_.grid->Ny; // dimension of linear system to solve
    size_t i, k;
    // when solving Zeta we fill linsys with dyy derivatives
    // Zeta.v is then solved exploiting normal Dirichlet boundary conditions
    // Zeta.u and Zeta.w are solved exploiting tangent Dirichlet boundary conditions

    LinearSys mySystem_u(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_v(sysDimension, BoundaryType::Normal);
    LinearSys mySystem_w(sysDimension, BoundaryType::Tangent);
    std::vector<double> rhs_u(sysDimension);
    std::vector<double> rhs_v(sysDimension);
    std::vector<double> rhs_w(sysDimension);

    std::vector<double> unknown_u;
    unknown_u.reserve(sysDimension);
    std::vector<double> unknown_v;
    unknown_v.reserve(sysDimension);
    std::vector<double> unknown_w;
    unknown_w.reserve(sysDimension);

    rhs_u.front() = 0;
    rhs_u.back() = 0;
    rhs_v.front() = 0;
    rhs_v.back() = 0;
    rhs_w.front() = 0;
    rhs_w.back() = 0;

    mul = 1.0 / (data_.grid->dy * data_.grid->dy);

    jStart = 0;
    for (i = 1; i < data_.grid->Nx-1; i++)
    {
        for (k = 1; k < data_.grid->Nz-1; k++)
        {
            iStart = i;
            kStart = k;

            for (size_t j = 1; j < sysDimension-1; j++)
            {
                porosity = data_.k(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                deriv_u = (data_.zeta(Axis::X, i, j + 1, k) + data_.zeta(Axis::X, i, j - 1, k) - 2.0 * data_.zeta(Axis::X, i, j, k))*mul;
                deriv_v = (data_.zeta(Axis::Y, i, j + 1, k) + data_.zeta(Axis::Y, i, j - 1, k) - 2.0 * data_.zeta(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.zeta(Axis::Z, i, j + 1, k) + data_.zeta(Axis::Z, i, j - 1, k) - 2.0 * data_.zeta(Axis::Z, i, j, k))*mul;

                rhs_u[j] = data_.eta(Axis::X, i,j,k) - gamma * deriv_u;
                rhs_v[j] = data_.eta(Axis::Y, i,j,k) - gamma * deriv_v;
                rhs_w[j] = data_.eta(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::X, Axis::Y, iStart, jStart, kStart);
            // mySystem_u.ThomaSolver();
            // std::vector<double> unknown_u = mySystem_u.getSolution();
            unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Y, Axis::Y, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v, BoundaryType::Normal);


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Z, Axis::Y, iStart, jStart, kStart);
            // mySystem_w.ThomaSolver();
            //std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);
            
            for (size_t j = 0; j < sysDimension; j++)
            {
                data_.zeta(Axis::X, i, j, k) = unknown_u[j];
                data_.zeta(Axis::Y, i, j, k) = unknown_v[j];
                data_.zeta(Axis::Z, i, j, k) = unknown_w[j];
            }

        }
    }

    // --- Handle Face Lines (Edges) ---
    // These loops run from 1 to N-2, avoiding the corner lines.
    // The logic is asymmetric (min vs max boundary).
    // Boundary values are always taken from data_.uBoundNew.

    // Special case, k = 0 (Z-min boundary)
    // Solve for Zeta.w (normal to Z), set Zeta.u, Zeta.v (tangents) from uBoundNew.
    k = 0;
    for (i = 1; i < data_.grid->Nx - 1; i++)
    {
        iStart = i;
        kStart = k;
        for (size_t j = 1; j < sysDimension - 1; j++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_w = (data_.zeta(Axis::Z, i, j + 1, k) + data_.zeta(Axis::Z, i, j - 1, k) - 2.0 * data_.zeta(Axis::Z, i, j, k)) * mul;
            // Note: rhs uses eta, as it's the input to this step
            rhs_w[j] = data_.eta(Axis::Z, i, j, k) - gamma * deriv_w;
        }

        mySystem_w.setRhs(rhs_w);
        mySystem_w.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Z, Axis::Y, iStart, jStart, kStart);
        //mySystem_w.ThomaSolver();
        //std::vector<double> unknown_w = mySystem_w.getSolution();
        unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);

        for (size_t j = 0; j < sysDimension; j++)
        {
            // Set tangents from boundary conditions
            data_.zeta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
            data_.zeta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
            // Set normal from solver
            data_.zeta(Axis::Z, i, j, k) = unknown_w[j];
        }
    }

    // Special case, k = Nz-1 (Z-max boundary)
    // Solve for Zeta.u, Zeta.v  (tangents to Z), set Zeta.w, (normal) from uBoundNew.
    k = data_.grid->Nz - 1;
    for (i = 1; i < data_.grid->Nx - 1; i++)
    {
        iStart = i;
        kStart = k;
        for (size_t j = 1; j < sysDimension - 1; j++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_u = (data_.zeta(Axis::X, i, j + 1, k) + data_.zeta(Axis::X, i, j - 1, k) - 2.0 * data_.zeta(Axis::X, i, j, k)) * mul;
            rhs_u[j] = data_.eta(Axis::X, i, j, k) - gamma * deriv_u;
            deriv_v = (data_.zeta(Axis::Y, i, j + 1, k) + data_.zeta(Axis::Y, i, j - 1, k) - 2.0 * data_.zeta(Axis::Y, i, j, k)) * mul;
            rhs_v[j] = data_.eta(Axis::Y, i, j, k) - gamma * deriv_v;
        }

        mySystem_u.setRhs(rhs_u);
        mySystem_u.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::X, Axis::Y, iStart, jStart, kStart);
        unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);

        mySystem_v.setRhs(rhs_v);
        mySystem_v.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Y, Axis::Y, iStart, jStart, kStart);
        unknown_v = solveSystem(mySystem_v, BoundaryType::Normal);


        for (size_t j = 0; j < sysDimension; j++)
        {
            data_.zeta(Axis::X, i, j, k) = unknown_u[j];
            data_.zeta(Axis::Y, i, j, k) = unknown_v[j];
            data_.zeta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
        }
    }

    // Special case, i = 0 (X-min boundary)
    // Solve for Zeta.u (normal to X), set Zeta.v, Zeta.w (tangents) from uBoundNew.
    i = 0;
    for (k = 1; k < data_.grid->Nz - 1; k++)
    {
        iStart = i;
        kStart = k;
        for (size_t j = 1; j < sysDimension - 1; j++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_u = (data_.zeta(Axis::X, i, j + 1, k) + data_.zeta(Axis::X, i, j - 1, k) - 2.0 * data_.zeta(Axis::X, i, j, k)) * mul;
            rhs_u[j] = data_.eta(Axis::X, i, j, k) - gamma * deriv_u;
            
        }

        mySystem_u.setRhs(rhs_u);
        mySystem_u.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::X, Axis::Y, iStart, jStart, kStart);
        //mySystem_u.ThomaSolver();
        //std::vector<double> unknown_u = mySystem_u.getSolution();
        unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);

        for (size_t j = 0; j < sysDimension; j++)
        {
            // Set normal from solver
            data_.zeta(Axis::X, i, j, k) = unknown_u[j];
            // Set tangents from boundary conditions
            data_.zeta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
            data_.zeta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
        }
    }

    // Special case, i = Nx-1 (X-max boundary)
    // Solve for Zeta.v (tangents to X), set Zeta.u,Zeta.w (normal) from uBoundNew.
    i = data_.grid->Nx - 1;
    for (k = 1; k < data_.grid->Nz - 1; k++)
    {
        iStart = i;
        kStart = k;
        for (size_t j = 1; j < sysDimension - 1; j++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_v = (data_.zeta(Axis::Y, i, j + 1, k) + data_.zeta(Axis::Y, i, j - 1, k) - 2.0 * data_.zeta(Axis::Y, i, j, k)) * mul;
            rhs_v[j] = data_.eta(Axis::Y, i, j, k) - gamma * deriv_v;
            deriv_w = (data_.zeta(Axis::Z, i, j + 1, k) + data_.zeta(Axis::Z, i, j - 1, k) - 2.0 * data_.zeta(Axis::Z, i, j, k)) * mul;
            rhs_w[j] = data_.eta(Axis::Z, i, j, k) - gamma * deriv_w;
        }
        mySystem_v.setRhs(rhs_v);
        mySystem_v.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Y, Axis::Y, iStart, jStart, kStart);
        unknown_v = solveSystem(mySystem_v, BoundaryType::Normal);

        mySystem_w.setRhs(rhs_w);
        mySystem_w.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Z, Axis::Y, iStart, jStart, kStart);
        unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);

        for (size_t j = 0; j < sysDimension; j++)
        {
            data_.zeta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
            data_.zeta(Axis::Y, i, j, k) = unknown_v[j];
            data_.zeta(Axis::Z, i, j, k) = unknown_w[j];
        }
    }


    // --- Handle Corner Lines (where i and k are boundaries) ---
    // For these lines, no solver is run. We apply boundary conditions directly.
    // Note: sysDimension == data_.grid->Ny here.

    // Special case: Corner i = 0, k = 0
    i = 0;
    k = 0;
    for (size_t j = 0; j < sysDimension; j++) // Loop is over j
    {
        data_.zeta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
        data_.zeta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
        data_.zeta(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
    }

    // Special case: Corner i = Nx-1, k = Nz-1
    i = data_.grid->Nx - 1;
    k = data_.grid->Nz - 1;
    for (size_t j = 1; j < sysDimension - 1; j++)
    {
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_v = (data_.zeta(Axis::Y, i, j + 1, k) + data_.zeta(Axis::Y, i, j - 1, k) - 2.0 * data_.zeta(Axis::Y, i, j, k)) * mul;
        rhs_v[j] = data_.eta(Axis::Y, i, j, k) - gamma * deriv_v;
    }
    mySystem_v.setRhs(rhs_v);
    mySystem_v.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Y, Axis::Y, iStart, jStart, kStart);
    unknown_v = solveSystem(mySystem_v, BoundaryType::Normal);

    for (size_t j = 0; j < sysDimension; j++) // Loop is over j
    {
        data_.zeta(Axis::X, i, j, k) = data_.bcu((i + 0.5) * dx, j * dy, k * dz, t);
        data_.zeta(Axis::Y, i, j, k) = unknown_v[j];
        data_.zeta(Axis::Z, i, j, k) = data_.bcw(i * dx, j * dy, (k + 0.5) * dz, t);
    }

    // Special case: Corner i = 0, k = Nz-1
    i = 0;
    k = data_.grid->Nz - 1;
    for (size_t j = 1; j < sysDimension - 1; j++)
    {
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_u = (data_.zeta(Axis::X, i, j + 1, k) + data_.zeta(Axis::X, i, j - 1, k) - 2.0 * data_.zeta(Axis::X, i, j, k)) * mul;
        rhs_u[j] = data_.eta(Axis::X, i, j, k) - gamma * deriv_u;
    }

    mySystem_u.setRhs(rhs_u);
    mySystem_u.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::X, Axis::Y, iStart, jStart, kStart);
    unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);

    for (size_t j = 0; j < sysDimension; j++) // Loop is over j
    {
        data_.zeta(Axis::X, i, j, k) = unknown_u[j];
        data_.zeta(Axis::Y, i, j, k) = data_.bcv(i * dx, (j + 0.5) * dy, k * dz, t);
        data_.zeta(Axis::Z, i, j, k) = data_.bcw(i * dx, j * dy, (k + 0.5) * dz, t);
    }

    // Special case: Corner i = Nx-1, k = 0
    i = data_.grid->Nx - 1;
    k = 0;
    for (size_t j = 1; j < sysDimension - 1; j++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_w = (data_.zeta(Axis::Z, i, j + 1, k) + data_.zeta(Axis::Z, i, j - 1, k) - 2.0 * data_.zeta(Axis::Z, i, j, k)) * mul;
            rhs_w[j] = data_.eta(Axis::Z, i, j, k) - gamma * deriv_w;
        }
        mySystem_w.setRhs(rhs_w);
        mySystem_w.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Z, Axis::Y, iStart, jStart, kStart);
        unknown_w = solveSystem(mySystem_w, BoundaryType::Tangent);

    for (size_t j = 0; j < sysDimension; j++) // Loop is over j
    {
        data_.zeta(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
        data_.zeta(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
        data_.zeta(Axis::Z, i, j, k) = unknown_w[j];
    }


    }





    // ------------------------------------------
    // Solve U --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::Z;
    size_t nSystem = data_.grid->Nx * data_.grid->Ny; // number of linear systems to solve
    size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve
    size_t i, j;
    // when solving U we fill linsys with dzz derivatives
    // U.w is then solved exploiting normal Dirichlet boundary conditions
    // U.v and U.u are solved exploiting tangent Dirichlet boundary conditions

    LinearSys mySystem_u(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_v(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_w(sysDimension, BoundaryType::Normal);
    std::vector<double> rhs_u(sysDimension);
    std::vector<double> rhs_v(sysDimension);
    std::vector<double> rhs_w(sysDimension);

    std::vector<double> unknown_u;
    unknown_u.reserve(sysDimension);
    std::vector<double> unknown_v;
    unknown_v.reserve(sysDimension);
    std::vector<double> unknown_w;
    unknown_w.reserve(sysDimension);

    rhs_u.front() = 0;
    rhs_u.back() = 0;
    rhs_v.front() = 0;
    rhs_v.back() = 0;
    rhs_w.front() = 0;
    rhs_w.back() = 0;

    mul = 1.0 / (data_.grid->dz * data_.grid->dz);
    kStart = 0;
    for (i = 1; i < data_.grid->Nx-1; i++)
    {
        for (j = 1; j < data_.grid->Ny-1; j++)
        {
            iStart = i;
            jStart = j;

            for (size_t k = 1; k < sysDimension-1; k++)
            {
                porosity = data_.k(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                deriv_u = (data_.u(Axis::X, i, j, k + 1) + data_.u(Axis::X, i, j, k - 1) - 2.0 * data_.u(Axis::X, i, j, k))*mul;
                deriv_v = (data_.u(Axis::Y, i, j, k + 1) + data_.u(Axis::Y, i, j, k - 1) - 2.0 * data_.u(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.u(Axis::Z, i, j, k + 1) + data_.u(Axis::Z, i, j, k - 1) - 2.0 * data_.u(Axis::Z, i, j, k))*mul;

                rhs_u[k] = data_.zeta(Axis::X, i,j,k) - gamma * deriv_u;
                rhs_v[k] = data_.zeta(Axis::Y, i,j,k) - gamma * deriv_v;
                rhs_w[k] = data_.zeta(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::X, Axis::Z, iStart, jStart, kStart);
            //mySystem_u.ThomaSolver();
            //std::vector<double> unknown_u = mySystem_u.getSolution();
            unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Y, Axis::Z, iStart, jStart, kStart);
            //mySystem_v.ThomaSolver();
            //std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Z, Axis::Z, iStart, jStart, kStart);
            //mySystem_w.ThomaSolver();
            //std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w, BoundaryType::Normal);
            
            for (size_t k = 0; k < sysDimension; k++)
            {
                data_.u(Axis::X, i, j, k) = unknown_u[k];
                data_.u(Axis::Y, i, j, k) = unknown_v[k];
                data_.u(Axis::Z, i, j, k) = unknown_w[k];
            }

        }
    }
    // --- Handle Face Lines (Edges) ---
    // These loops run from 1 to N-2, avoiding the corner lines.
    // The logic is asymmetric (min vs max boundary).
    // Boundary values are always taken from data_.uBoundNew.

    // Special case, j = 0 (Y-min boundary)
    // Solve for U.v (normal to Y), set U.u, U.w (tangents) from uBoundNew.
    j = 0;
    for (i = 1; i < data_.grid->Nx - 1; i++)
    {
        iStart = i;
        jStart = j;
        for (size_t k = 1; k < sysDimension - 1; k++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_v = (data_.u(Axis::Y, i, j, k + 1) + data_.u(Axis::Y, i, j, k - 1) - 2.0 * data_.u(Axis::Y, i, j, k)) * mul;
            // Note: rhs uses zeta, as it's the input to this step
            rhs_v[k] = data_.zeta(Axis::Y, i, j, k) - gamma * deriv_v;
        }

        mySystem_v.setRhs(rhs_v);
        mySystem_v.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Y, Axis::Z, iStart, jStart, kStart);
        //mySystem_v.ThomaSolver();
        //std::vector<double> unknown_v = mySystem_v.getSolution();
        unknown_v = solveSystem(mySystem_v, BoundaryType::Normal);

        for (size_t k = 0; k < sysDimension; k++)
        {
            data_.u(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
            data_.u(Axis::Y, i, j, k) = unknown_v[k];
            data_.u(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
        }
    }

    // Special case, j = Ny-1 (Y-max boundary)
    // Solve for U.u (tangents to Y), set U.v, U.w (normal) from uBoundNew.
    j = data_.grid->Ny - 1;
    for (i = 1; i < data_.grid->Nx - 1; i++)
    {
        iStart = i;
        jStart = j;
        for (size_t k = 1; k < sysDimension - 1; k++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_u = (data_.u(Axis::X, i, j, k + 1) + data_.u(Axis::X, i, j, k - 1) - 2.0 * data_.u(Axis::X, i, j, k)) * mul;
            rhs_u[k] = data_.zeta(Axis::X, i, j, k) - gamma * deriv_u;
            deriv_w = (data_.u(Axis::Z, i, j, k + 1) + data_.u(Axis::Z, i, j, k - 1) - 2.0 * data_.u(Axis::Z, i, j, k)) * mul;
            rhs_w[k] = data_.zeta(Axis::Z, i, j, k) - gamma * deriv_w;
        }

        mySystem_u.setRhs(rhs_u);
        mySystem_u.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::X, Axis::Z, iStart, jStart, kStart);
        unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);

        mySystem_w.setRhs(rhs_w);
        mySystem_w.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Z, Axis::Z, iStart, jStart, kStart);
        unknown_w = solveSystem(mySystem_w, BoundaryType::Normal);

        for (size_t k = 0; k < sysDimension; k++)
        {
            // Set tangents from solver
            data_.u(Axis::X, i, j, k) = unknown_u[k];
            data_.u(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
            data_.u(Axis::Z, i, j, k) = unknown_w[k];
            
        }
    }

    // Special case, i = 0 (X-min boundary)
    // Solve for U.u (normal to X), set U.v, U.w (tangents) from uBoundNew.
    i = 0;
    for (j = 1; j < data_.grid->Ny - 1; j++)
    {
        iStart = i;
        jStart = j;
        for (size_t k = 1; k < sysDimension - 1; k++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_u = (data_.u(Axis::X, i, j, k + 1) + data_.u(Axis::X, i, j, k - 1) - 2.0 * data_.u(Axis::X, i, j, k)) * mul;
            rhs_u[k] = data_.zeta(Axis::X, i, j, k) - gamma * deriv_u;
        }

        mySystem_u.setRhs(rhs_u);
        mySystem_u.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::X, Axis::Z, iStart, jStart, kStart);
        //mySystem_u.ThomaSolver();
        //std::vector<double> unknown_u = mySystem_u.getSolution();
        unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);

        for (size_t k = 0; k < sysDimension; k++)
        {
            // Set normal from solver
            data_.u(Axis::X, i, j, k) = unknown_u[k];
            data_.u(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
            data_.u(Axis::Z, i, j, k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, t);
        }
    }

    // Special case, i = Nx-1 (X-max boundary)
    // Solve for U.v (tangents to X), set U.u, U.w  (normal) from uBoundNew.
    i = data_.grid->Nx - 1;
    for (j = 1; j < data_.grid->Ny - 1; j++)
    {
        iStart = i;
        jStart = j;
        for (size_t k = 1; k < sysDimension - 1; k++)
        {
            porosity = data_.k(i, j, k);
            beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
            gamma = data_.dt * data_.nu * 0.5 / beta;
            deriv_v = (data_.u(Axis::Y, i, j, k + 1) + data_.u(Axis::Y, i, j, k - 1) - 2.0 * data_.u(Axis::Y, i, j, k)) * mul;
            rhs_v[k] = data_.zeta(Axis::Y, i, j, k) - gamma * deriv_v;
            deriv_w = (data_.u(Axis::Z, i, j, k + 1) + data_.u(Axis::Z, i, j, k - 1) - 2.0 * data_.u(Axis::Z, i, j, k)) * mul;
            rhs_w[k] = data_.zeta(Axis::Z, i, j, k) - gamma * deriv_w;
        }

        mySystem_v.setRhs(rhs_v);
        mySystem_v.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Y, Axis::Z, iStart, jStart, kStart);
        unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);
        mySystem_w.setRhs(rhs_w);
        mySystem_w.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Z, Axis::Z, iStart, jStart, kStart);
        unknown_w = solveSystem(mySystem_w, BoundaryType::Normal);

        for (size_t k = 0; k < sysDimension; k++)
        {
            // Set normal from boundary conditions
            data_.u(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
            data_.u(Axis::Y, i, j, k) = unknown_v[k];
            data_.u(Axis::Z, i, j, k) = unknown_w[k];
        }
    }


    // --- Handle Corner Lines (where i and j are boundaries) ---
    // For these lines, no solver is run. We apply boundary conditions directly.
    // Note: sysDimension == data_.grid->Nz here.

    // Special case: Corner i = 0, j = 0
    i = 0;
    j = 0;
    for (size_t k = 0; k < sysDimension; k++) // Loop is over k
    {
        data_.u(Axis::X, i, j, k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, t);
        data_.u(Axis::Y, i, j, k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, t);
        data_.u(Axis::Z, i, j, k) = data_.bcw(i * dx, j * dy, (k + 0.5) * dz, t);
    }

    // Special case: Corner i = Nx-1, j = Ny-1
    i = data_.grid->Nx - 1;
    j = data_.grid->Ny - 1;
    for (size_t k = 1; k < sysDimension - 1; k++)
    {
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_w = (data_.u(Axis::Z, i, j, k + 1) + data_.u(Axis::Z, i, j, k - 1) - 2.0 * data_.u(Axis::Z, i, j, k)) * mul;
        rhs_w[k] = data_.zeta(Axis::Z, i, j, k) - gamma * deriv_w;
    }

    mySystem_w.setRhs(rhs_w);
    mySystem_w.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Z, Axis::Z, iStart, jStart, kStart);
    unknown_w = solveSystem(mySystem_w, BoundaryType::Normal);

    for (size_t k = 0; k < sysDimension; k++) // Loop is over k
    {
        data_.u(Axis::X, i, j, k) = data_.bcu((i + 0.5) * dx, j * dy, k * dz, t);
        data_.u(Axis::Y, i, j, k) = data_.bcv(i * dx, (j + 0.5) * dy, k * dz, t);
        data_.u(Axis::Z, i, j, k) = unknown_w[k];
    }

    // Special case: Corner i = 0, j = Ny-1
    i = 0;
    j = data_.grid->Ny - 1;
    for (size_t k = 1; k < sysDimension - 1; k++)
    {
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_u = (data_.u(Axis::X, i, j, k + 1) + data_.u(Axis::X, i, j, k - 1) - 2.0 * data_.u(Axis::X, i, j, k)) * mul;
        rhs_u[k] = data_.zeta(Axis::X, i, j, k) - gamma * deriv_u;
    }

    mySystem_u.setRhs(rhs_u);
    mySystem_u.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::X, Axis::Z, iStart, jStart, kStart);
    unknown_u = solveSystem(mySystem_u, BoundaryType::Tangent);

    for (size_t k = 0; k < sysDimension; k++) // Loop is over k
    {
        data_.u(Axis::X, i, j, k) = unknown_u[k];
        data_.u(Axis::Y, i, j, k) = data_.bcv(i * dx, (j + 0.5) * dy, k * dz, t);
        data_.u(Axis::Z, i, j, k) = data_.bcw(i * dx, j * dy, (k + 0.5) * dz, t);
    }

    // Special case: Corner i = Nx-1, j = 0
    i = data_.grid->Nx - 1;
    j = 0;
    for (size_t k = 1; k < sysDimension - 1; k++)
    {
        porosity = data_.k(i, j, k);
        beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
        gamma = data_.dt * data_.nu * 0.5 / beta;
        deriv_v = (data_.u(Axis::Y, i, j, k + 1) + data_.u(Axis::Y, i, j, k - 1) - 2.0 * data_.u(Axis::Y, i, j, k)) * mul;
        rhs_v[k] = data_.zeta(Axis::Y, i, j, k) - gamma * deriv_v;
    }

    mySystem_v.setRhs(rhs_v);
    mySystem_v.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Y, Axis::Z, iStart, jStart, kStart);
    unknown_v = solveSystem(mySystem_v, BoundaryType::Tangent);

    for (size_t k = 0; k < sysDimension; k++) // Loop is over k
    {
        data_.u(Axis::X, i, j, k) = data_.bcu((i + 0.5) * dx, j * dy, k * dz, t);
        data_.u(Axis::Y, i, j, k) = unknown_v[k];
        data_.u(Axis::Z, i, j, k) = data_.bcw(i * dx, j * dy, (k + 0.5) * dz, t);
    }
    }
}
