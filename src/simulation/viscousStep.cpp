#include <algorithm>
#include <cmath>
#include <core/Fields.hpp>
#include <vector>
#include <numerics/derivatives.hpp>
#include <numerics/LinearSys.hpp>
#include "numerics/SchurSequentialSolver.hpp"
#include <simulation/viscousStep.hpp>


std::vector<double> ViscousStep::solveSystem(LinearSys& sys) {

    // 1. If we want to use the classic solver (debug or actual P=1)
    if (parallel_.schurDomains <= 1) {
        sys.ThomaSolver();
        return sys.getSolution();
    }
    else
    {
        // 2. Sequential Schur Logic
        SchurSequentialSolver schur(sys.getMatrix().getSize(), parallel_.schurDomains);
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
    auto &predictor = data_.predictor;
    double nu_val = data_.nu;

    Derivatives derive;
    derive.computeGradient(predictor, gradP);

    // Recepie
    // g = f  -grad(p)  -nu/k *u * 0.5  +nu*(dxx eta + dyy zeta + dzz u) * 0.5

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z})
    {

        for (size_t k = 0; k < data_.grid->Nz; k++)
            for (size_t j = 0; j < data_.grid->Ny; j++)
                for (size_t i = 0; i < data_.grid->Nx; i++)
                {
                    dxxEta(axis, i, j, k) = derive.Dxx_local(eta(axis), i, j, k);
                    dyyZeta(axis, i, j, k) = derive.Dyy_local(zeta(axis), i, j, k);
                    dzzU(axis, i, j, k) = derive.Dzz_local(u(axis), i, j, k);
                    
                }

        auto &f_data = data_.f(axis).getData();
        auto &u_data = data_.u(axis).getData();
        auto &gradP_data = gradP(axis).getData();
        auto &dxx_data = dxxEta(axis).getData();
        auto &dyy_data = dyyZeta(axis).getData();
        auto &dzz_data = dzzU(axis).getData();
        auto &g_data = g(axis).getData();
        auto &k_data = data_.k(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++)
        {
            g_data[i] = f_data[i] - gradP_data[i] - nu_val * u_data[i] / k_data[i] + nu_val * (dxx_data[i] + dyy_data[i] + dzz_data[i]);
        }
    }
    }

void ViscousStep::computeXi()
{
    // Ingredients list
    auto& u = data_.u;
    double nu_val = data_.nu;
    double dt_val = data_.dt;
   
    // Recepie
    // beta = 1+ dt*nu /2/k
    // xi = u + dt/beta * g

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {

        auto& u_data = data_.u(axis).getData();
        auto& k_data = data_.k(axis).getData();
        auto& g_data = g(axis).getData();
        auto& xi_data = xi(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++)
        {
            xi_data[i] = u_data[i]
                        + dt_val * g_data[i]
                        / (1 + dt_val * nu_val * 0.5 / k_data[i]);
        }
    }
    size_t i, j, k;
    auto &nx = data_.grid->Nx;
    auto &ny = data_.grid->Ny;
    auto &nz = data_.grid->Nz;
    auto &dx = data_.grid->dx;
    auto &dy = data_.grid->dy;
    auto &dz = data_.grid->dz;
    auto &time = data_.currTime;
    
    i = 0;
    for (j = 0; j < ny; j++)
    {
        for (k = 0; k < nz; k++)
        {
            xi(Axis::Y)(i,j,k) = data_.bcv(0.0, (j+0.5)*dy, k*dz, time);
            xi(Axis::Z)(i,j,k) = data_.bcw(0.0, j*dy, (k+0.5)*dz, time);
        }
    }
    j = 0;
    for (i = 0; i < nx; i++)
    {
        for (k = 0; k < nz; k++)
        {
            xi(Axis::X)(i,j,k) = data_.bcu((i+0.5)*dx, 0.0, k*dz, time);
            xi(Axis::Z)(i,j,k) = data_.bcw(i*dx, 0.0, (k+0.5)*dz, time);
        }
    }
    k = 0;
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            xi(Axis::Y)(i,j,k) = data_.bcv(i*dx, (j+0.5)*dy, 0.0, time);
            xi(Axis::X)(i,j,k) = data_.bcu((i+0.5)*dx, j*dy, 0.0, time);
        }
    }

    i = nx - 1;
    for (j = 0; j < ny; j++)
    {
        for (k = 0; k < nz; k++)
        {
            xi(Axis::X)(i,j,k) = data_.bcu((i+0.5)*dx, j*dy, k*dz, time);
        }
    }
    j = ny - 1;
    for (i = 0; i < nx; i++)
    {
        for (k = 0; k < nz; k++)
        {
            xi(Axis::Y)(i,j,k) = data_.bcv(i*dx, (j+0.5)*dy, k*dz, time);
        }
    }
    k = nz - 1;
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            xi(Axis::Z)(i,j,k) = data_.bcw(i*dx, j*dy, (k+0.5)*dz, time);
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

    LinearSys mySystem_u(sysDimension);
    LinearSys mySystem_v(sysDimension);
    LinearSys mySystem_w(sysDimension);
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

    for (j = 0; j < data_.grid->Ny; j++)
    {
        for (k = 0; k < data_.grid->Nz; k++)
        {
            jStart = j;
            kStart = k;

            for (size_t i = 1; i < sysDimension-1; i++)
            {
                deriv_u = (data_.eta(Axis::X, i + 1, j, k) + data_.eta(Axis::X, i - 1, j, k) - 2.0 * data_.eta(Axis::X, i, j, k))*mul;
                deriv_v = (data_.eta(Axis::Y, i + 1, j, k) + data_.eta(Axis::Y, i - 1, j, k) - 2.0 * data_.eta(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.eta(Axis::Z, i + 1, j, k) + data_.eta(Axis::Z, i - 1, j, k) - 2.0 * data_.eta(Axis::Z, i, j, k))*mul;


                porosity = data_.k(Axis::X)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_u[i] = xi(Axis::X, i,j,k) - gamma * deriv_u;

                porosity = data_.k(Axis::Y)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_v[i] = xi(Axis::Y, i,j,k) - gamma * deriv_v;

                porosity = data_.k(Axis::Z)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_w[i] = xi(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);
            
            mySystem_u.fillSystemVelocity(data_, data_.eta, xi, Axis::X, Axis::X, iStart, jStart, kStart);
            // mySystem_u.ThomaSolver();
            // std::vector<double> unknown_u = mySystem_u.getSolution();
            unknown_u = solveSystem(mySystem_u);


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.eta, xi, Axis::Y, Axis::X, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v);


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.eta, xi, Axis::Z, Axis::X, iStart, jStart, kStart);
            // mySystem_w.ThomaSolver();
            // std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w);
            
            for (size_t i = 0; i < sysDimension; i++)
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
    normalAxis = Axis::Y;
    size_t nSystem = data_.grid->Nx * data_.grid->Nz; // number of linear systems to solve
    size_t sysDimension = data_.grid->Ny; // dimension of linear system to solve
    size_t i, k;
    // when solving Zeta we fill linsys with dyy derivatives
    // Zeta.v is then solved exploiting normal Dirichlet boundary conditions
    // Zeta.u and Zeta.w are solved exploiting tangent Dirichlet boundary conditions

    LinearSys mySystem_u(sysDimension);
    LinearSys mySystem_v(sysDimension);
    LinearSys mySystem_w(sysDimension);
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
    for (i = 0; i < data_.grid->Nx; i++)
    {
        for (k = 0; k < data_.grid->Nz; k++)
        {
            iStart = i;
            kStart = k;

            for (size_t j = 1; j < sysDimension-1; j++)
            {
                deriv_u = (data_.zeta(Axis::X, i, j + 1, k) + data_.zeta(Axis::X, i, j - 1, k) - 2.0 * data_.zeta(Axis::X, i, j, k))*mul;
                deriv_v = (data_.zeta(Axis::Y, i, j + 1, k) + data_.zeta(Axis::Y, i, j - 1, k) - 2.0 * data_.zeta(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.zeta(Axis::Z, i, j + 1, k) + data_.zeta(Axis::Z, i, j - 1, k) - 2.0 * data_.zeta(Axis::Z, i, j, k))*mul;

                porosity = data_.k(Axis::X)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_u[j] = data_.eta(Axis::X, i,j,k) - gamma * deriv_u;

                porosity = data_.k(Axis::Y)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_v[j] = data_.eta(Axis::Y, i,j,k) - gamma * deriv_v;

                porosity = data_.k(Axis::Z)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_w[j] = data_.eta(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::X, Axis::Y, iStart, jStart, kStart);
            // mySystem_u.ThomaSolver();
            // std::vector<double> unknown_u = mySystem_u.getSolution();
            unknown_u = solveSystem(mySystem_u);


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Y, Axis::Y, iStart, jStart, kStart);
            // mySystem_v.ThomaSolver();
            // std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v);


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.zeta, data_.eta, Axis::Z, Axis::Y, iStart, jStart, kStart);
            // mySystem_w.ThomaSolver();
            //std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w);
            
            for (size_t j = 0; j < sysDimension; j++)
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
    normalAxis = Axis::Z;
    size_t nSystem = data_.grid->Nx * data_.grid->Ny; // number of linear systems to solve
    size_t sysDimension = data_.grid->Nz; // dimension of linear system to solve
    size_t i, j;
    // when solving U we fill linsys with dzz derivatives
    // U.w is then solved exploiting normal Dirichlet boundary conditions
    // U.v and U.u are solved exploiting tangent Dirichlet boundary conditions

    LinearSys mySystem_u(sysDimension);
    LinearSys mySystem_v(sysDimension);
    LinearSys mySystem_w(sysDimension);
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
    for (i = 0; i < data_.grid->Nx; i++)
    {
        for (j = 0; j < data_.grid->Ny; j++)
        {
            iStart = i;
            jStart = j;

            for (size_t k = 1; k < sysDimension-1; k++)
            {
                deriv_u = (data_.u(Axis::X, i, j, k + 1) + data_.u(Axis::X, i, j, k - 1) - 2.0 * data_.u(Axis::X, i, j, k))*mul;
                deriv_v = (data_.u(Axis::Y, i, j, k + 1) + data_.u(Axis::Y, i, j, k - 1) - 2.0 * data_.u(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.u(Axis::Z, i, j, k + 1) + data_.u(Axis::Z, i, j, k - 1) - 2.0 * data_.u(Axis::Z, i, j, k))*mul;

                porosity = data_.k(Axis::X)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_u[k] = data_.zeta(Axis::X, i,j,k) - gamma * deriv_u;

                porosity = data_.k(Axis::Y)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_v[k] = data_.zeta(Axis::Y, i,j,k) - gamma * deriv_v;

                porosity = data_.k(Axis::Z)(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; 
                rhs_w[k] = data_.zeta(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::X, Axis::Z, iStart, jStart, kStart);
            //mySystem_u.ThomaSolver();
            //std::vector<double> unknown_u = mySystem_u.getSolution();
            unknown_u = solveSystem(mySystem_u);


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Y, Axis::Z, iStart, jStart, kStart);
            //mySystem_v.ThomaSolver();
            //std::vector<double> unknown_v = mySystem_v.getSolution();
            unknown_v = solveSystem(mySystem_v);


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_, data_.u, data_.zeta, Axis::Z, Axis::Z, iStart, jStart, kStart);
            //mySystem_w.ThomaSolver();
            //std::vector<double> unknown_w = mySystem_w.getSolution();
            unknown_w = solveSystem(mySystem_w);
            
            for (size_t k = 0; k < sysDimension; k++)
            {
                data_.u(Axis::X, i, j, k) = unknown_u[k];
                data_.u(Axis::Y, i, j, k) = unknown_v[k];
                data_.u(Axis::Z, i, j, k) = unknown_w[k];
            }

        }
    }
    }
}
