#include <algorithm>
#include <cmath>
#include <core/Fields.hpp>
#include <vector>
#include <numerics/derivatives.hpp>
#include <numerics/LinearSys.hpp>
#include <simulation/viscousStep.hpp>



ViscousStep::ViscousStep(SimulationData& simData) : data_(simData)
{
    initializeWorkspaceFields(data_.gridPtr);
}



void ViscousStep::run()
{
    computeG();
    computeXi();
    closeViscousStep();
}

void ViscousStep::initializeWorkspaceFields(std::shared_ptr<const Grid> gridPtr)
{
    std::vector<Field::Scalar> zeros(gridPtr->size(), 0.0);
    g.setup(gridPtr, zeros, zeros, zeros);
    gradP.setup(gridPtr, zeros, zeros, zeros);
    dxxEta.setup(gridPtr, zeros, zeros, zeros);
    dyyZeta.setup(gridPtr, zeros, zeros, zeros);
    dzzU.setup(gridPtr, zeros, zeros, zeros);
    xi.setup(gridPtr, zeros, zeros, zeros);
}

void ViscousStep::computeG()
{
    // Ingredients list
    auto& eta = data_.eta;
    auto& zeta = data_.zeta;
    auto& u = data_.u;
    auto& p = data_.p;
    auto& k_data = data_.k.getData(); // k cant be 0!!!
    double nu_val = data_.nu;
    
    Derivatives derive;
    derive.computeGradient(p, gradP);
    

    // Recepie
    // g = f  -grad(p)  -nu/k *u * 0.5  +nu*(dxx eta + dyy zeta + dzz u) * 0.5

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {

        derive.computeDxx(eta(axis), dxxEta(axis));
        derive.computeDyy(zeta(axis), dyyZeta(axis));
        derive.computeDzz(u(axis), dzzU(axis));

        auto& f_data = data_.f(axis).getData();
        auto& u_data = data_.u(axis).getData();
        auto& gradP_data = gradP(axis).getData();
        auto& dxx_data = dxxEta(axis).getData();
        auto& dyy_data = dyyZeta(axis).getData();
        auto& dzz_data = dzzU(axis).getData();
        auto& g_data = g(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++)
        {
            g_data[i] = f_data[i] 
                        - gradP_data[i] 
                        - nu_val * u_data[i] / k_data[i] * 0.5
                        + nu_val * (dxx_data[i] + dyy_data[i] + dzz_data[i]) * 0.5;
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
    double porosity, beta, gamma, mul, deriv_u, deriv_v, deriv_w;
    size_t iStart, jStart, kStart;
    // solve on the face orthogonal to normalAxis
    Axis normalAxis;
    

    // ------------------------------------------
    // Solve Eta --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::X;
    size_t nSystem = data_.gridPtr->Ny * data_.gridPtr->Nz; // number of linear systems to solve
    size_t sysDimension = data_.gridPtr->Nx; // dimension of linear system to solve
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
    rhs_u.front() = 0;
    rhs_u.back() = 0;
    rhs_v.front() = 0;
    rhs_v.back() = 0;
    rhs_w.front() = 0;
    rhs_w.back() = 0;

    iStart = 0;
    for (size_t j = 1; j < data_.gridPtr->Ny-1; j++)
    {
        for (size_t k = 1; k < data_.gridPtr->Nz-1; k++)
        {
            jStart = j;
            kStart = k;

            for (size_t i = 1; i < sysDimension-1; i++)
            {
                porosity = data_.k(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                mul = 1.0 / (data_.gridPtr->dx * data_.gridPtr->dx);
                deriv_u = (data_.eta(Axis::X, i + 1, j, k) + data_.eta(Axis::X, i - 1, j, k) - 2.0 * data_.eta(Axis::X, i, j, k))*mul;
                deriv_v = (data_.eta(Axis::Y, i + 1, j, k) + data_.eta(Axis::Y, i - 1, j, k) - 2.0 * data_.eta(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.eta(Axis::Z, i + 1, j, k) + data_.eta(Axis::Z, i - 1, j, k) - 2.0 * data_.eta(Axis::Z, i, j, k))*mul;

                rhs_u[i] = xi(Axis::X, i,j,k) - gamma * deriv_u;
                rhs_v[i] = xi(Axis::Y, i,j,k) - gamma * deriv_v;
                rhs_w[i] = xi(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);
            
            mySystem_u.fillSystemVelocity(data_.k, data_.eta, xi, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::X, Axis::X, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_u.ThomaSolver();
            std::vector<double> unknown_u = mySystem_u.getSolution();


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_.k, data_.eta, xi, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::Y, Axis::X, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_v.ThomaSolver();
            std::vector<double> unknown_v = mySystem_v.getSolution();


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_.k, data_.eta, xi, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::Z, Axis::X, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_w.ThomaSolver();
            std::vector<double> unknown_w = mySystem_w.getSolution();
            
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
    size_t nSystem = data_.gridPtr->Nx * data_.gridPtr->Nz; // number of linear systems to solve
    size_t sysDimension = data_.gridPtr->Ny; // dimension of linear system to solve
    // when solving Zeta we fill linsys with dyy derivatives
    // Zeta.v is then solved exploiting normal Dirichlet boundary conditions
    // Zeta.u and Zeta.w are solved exploiting tangent Dirichlet boundary conditions

    LinearSys mySystem_u(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_v(sysDimension, BoundaryType::Normal);
    LinearSys mySystem_w(sysDimension, BoundaryType::Tangent);
    std::vector<double> rhs_u(sysDimension);
    std::vector<double> rhs_v(sysDimension);
    std::vector<double> rhs_w(sysDimension);
    rhs_u.front() = 0;
    rhs_u.back() = 0;
    rhs_v.front() = 0;
    rhs_v.back() = 0;
    rhs_w.front() = 0;
    rhs_w.back() = 0;


    jStart = 0;
    for (size_t i = 1; i < data_.gridPtr->Nx-1; i++)
    {
        for (size_t k = 1; k < data_.gridPtr->Nz-1; k++)
        {
            iStart = i;
            kStart = k;

            for (size_t j = 1; j < sysDimension-1; j++)
            {
                porosity = data_.k(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                mul = 1.0 / (data_.gridPtr->dy * data_.gridPtr->dy);
                deriv_u = (data_.zeta(Axis::X, i, j + 1, k) + data_.zeta(Axis::X, i, j - 1, k) - 2.0 * data_.zeta(Axis::X, i, j, k))*mul;
                deriv_v = (data_.zeta(Axis::Y, i, j + 1, k) + data_.zeta(Axis::Y, i, j - 1, k) - 2.0 * data_.zeta(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.zeta(Axis::Z, i, j + 1, k) + data_.zeta(Axis::Z, i, j - 1, k) - 2.0 * data_.zeta(Axis::Z, i, j, k))*mul;

                rhs_u[j] = data_.eta(Axis::X, i,j,k) - gamma * deriv_u;
                rhs_v[j] = data_.eta(Axis::Y, i,j,k) - gamma * deriv_v;
                rhs_w[j] = data_.eta(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_.k, data_.zeta, data_.eta, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::X, Axis::Y, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_u.ThomaSolver();
            std::vector<double> unknown_u = mySystem_u.getSolution();


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_.k, data_.zeta, data_.eta, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::Y, Axis::Y, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_v.ThomaSolver();
            std::vector<double> unknown_v = mySystem_v.getSolution();


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_.k, data_.zeta, data_.eta, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::Z, Axis::Y, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_w.ThomaSolver();
            std::vector<double> unknown_w = mySystem_w.getSolution();
            
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
    size_t nSystem = data_.gridPtr->Nx * data_.gridPtr->Ny; // number of linear systems to solve
    size_t sysDimension = data_.gridPtr->Nz; // dimension of linear system to solve
    // when solving U we fill linsys with dzz derivatives
    // U.w is then solved exploiting normal Dirichlet boundary conditions
    // U.v and U.u are solved exploiting tangent Dirichlet boundary conditions

    LinearSys mySystem_u(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_v(sysDimension, BoundaryType::Tangent);
    LinearSys mySystem_w(sysDimension, BoundaryType::Normal);
    std::vector<double> rhs_u(sysDimension);
    std::vector<double> rhs_v(sysDimension);
    std::vector<double> rhs_w(sysDimension);
    rhs_u.front() = 0;
    rhs_u.back() = 0;
    rhs_v.front() = 0;
    rhs_v.back() = 0;
    rhs_w.front() = 0;
    rhs_w.back() = 0;


    kStart = 0;
    for (size_t i = 1; i < data_.gridPtr->Nx-1; i++)
    {
        for (size_t j = 1; j < data_.gridPtr->Ny-1; j++)
        {
            iStart = i;
            jStart = j;

            for (size_t k = 1; k < sysDimension-1; k++)
            {
                porosity = data_.k(i,j,k);
                beta = 1 + (data_.dt * data_.nu * 0.5 / porosity);
                gamma = data_.dt * data_.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                mul = 1.0 / (data_.gridPtr->dz * data_.gridPtr->dz);
                deriv_u = (data_.u(Axis::X, i, j, k + 1) + data_.u(Axis::X, i, j, k - 1) - 2.0 * data_.u(Axis::X, i, j, k))*mul;
                deriv_v = (data_.u(Axis::Y, i, j, k + 1) + data_.u(Axis::Y, i, j, k - 1) - 2.0 * data_.u(Axis::Y, i, j, k))*mul;
                deriv_w = (data_.u(Axis::Z, i, j, k + 1) + data_.u(Axis::Z, i, j, k - 1) - 2.0 * data_.u(Axis::Z, i, j, k))*mul;

                rhs_u[k] = data_.zeta(Axis::X, i,j,k) - gamma * deriv_u;
                rhs_v[k] = data_.zeta(Axis::Y, i,j,k) - gamma * deriv_v;
                rhs_w[k] = data_.zeta(Axis::Z, i,j,k) - gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(data_.k, data_.u, data_.zeta, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::X, Axis::Z, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_u.ThomaSolver();
            std::vector<double> unknown_u = mySystem_u.getSolution();


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(data_.k, data_.u, data_.zeta, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::Y, Axis::Z, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_v.ThomaSolver();
            std::vector<double> unknown_v = mySystem_v.getSolution();


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(data_.k, data_.u, data_.zeta, data_.uBoundNew, 
                                        data_.uBoundOld, Axis::Z, Axis::Z, iStart, jStart, kStart, data_.nu, data_.dt);
            mySystem_w.ThomaSolver();
            std::vector<double> unknown_w = mySystem_w.getSolution();
            
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







