#include <algorithm>
#include <cmath>
#include <simulation/viscousStep.hpp>
#include <core/Fields.hpp>
#include <vector>
#include <numerics/derivatives.hpp>
#include <numerics/LinearSys.hpp>



ViscousStep::ViscousStep(SimulationContext& ctx) : context_(ctx)
{
    initializeWorkspaceFields(context_.gridPtr);
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
    auto& eta = context_.state.etaOld;
    auto& zeta = context_.state.zetaOld;
    auto& u = context_.state.uOld;
    auto& p = context_.state.p;
    auto& k_data = context_.constants.k.getData(); // k cant be 0!!!
    double nu_val = context_.constants.nu;
    
    Derivatives derive;
    derive.computeGradient(p, gradP);
    

    // Recepie
    // g = f  -grad(p)  -nu/k *u  +nu*(dxx eta + dyy zeta + dzz u)

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {

        derive.computeDxx(eta(axis), dxxEta(axis));
        derive.computeDyy(zeta(axis), dyyZeta(axis));
        derive.computeDzz(u(axis), dzzU(axis));

        auto& f_data = context_.constants.f(axis).getData();
        auto& u_data = context_.state.uOld(axis).getData();
        auto& gradP_data = gradP(axis).getData();
        auto& dxx_data = dxxEta(axis).getData();
        auto& dyy_data = dyyZeta(axis).getData();
        auto& dzz_data = dzzU(axis).getData();
        auto& g_data = g(axis).getData();

        for (size_t i = 0; i < u_data.size(); i++)
        {
            g_data[i] = f_data[i] 
                        - gradP_data[i] 
                        - nu_val * u_data[i] / k_data[i]
                        + nu_val * (dxx_data[i] + dyy_data[i] + dzz_data[i]);
        }
    }
}

void ViscousStep::computeXi()
{
     // Ingredients list
    auto& u = context_.state.uOld;
    auto& k_data = context_.constants.k.getData();  // k cant be 0!!!
    double nu_val = context_.constants.nu;
    double dt_val = context_.timeSettings.dt;
   
    // Recepie
    // beta = 1+ dt*nu /2/k
    // xi = u + dt/beta * g

    // Let me cook
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {

        auto& u_data = context_.state.uOld(axis).getData();
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
    size_t nSystem = context_.gridPtr->Ny * context_.gridPtr->Nz; // number of linear systems to solve
    size_t sysDimension = context_.gridPtr->Nx; // dimension of linear system to solve
    // when solving Eta we fill linsys with dxx derivatives
    // Eta.u is then solved exploiting normal Neumann boundary conditions
    // Eta.v and Eta.w are solved exploiting tangent Dirichlet boundary conditions

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
    for (size_t j = 0; j < context_.gridPtr->Ny; j++)
    {
        for (size_t k = 0; k < context_.gridPtr->Nz; k++)
        {
            jStart = j;
            kStart = k;

            for (size_t i = 1; i < sysDimension-1; i++)
            {
                porosity = context_.constants.k(i,j,k);
                beta = 1 + (context_.timeSettings.dt * context_.constants.nu * 0.5 / porosity);
                gamma = context_.timeSettings.dt * context_.constants.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                mul = 1.0 / (context_.gridPtr->dx * context_.gridPtr->dx);
                deriv_u = (context_.state.etaOld(Axis::X, i + 1, j, k) + context_.state.etaOld(Axis::X, i - 1, j, k) - 2.0 * context_.state.etaOld(Axis::X, i, j, k))*mul;
                deriv_v = (context_.state.etaOld(Axis::Y, i + 1, j, k) + context_.state.etaOld(Axis::Y, i - 1, j, k) - 2.0 * context_.state.etaOld(Axis::Y, i, j, k))*mul;
                deriv_w = (context_.state.etaOld(Axis::Z, i + 1, j, k) + context_.state.etaOld(Axis::Z, i - 1, j, k) - 2.0 * context_.state.etaOld(Axis::Z, i, j, k))*mul;

                rhs_u[i] = xi(Axis::X, i,j,k) + gamma * deriv_u;
                rhs_v[i] = xi(Axis::Y, i,j,k) + gamma * deriv_v;
                rhs_w[i] = xi(Axis::Z, i,j,k) + gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(context_.constants.k, context_.state.etaOld, context_.state.xi, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::X, Axis::X, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_u.ThomaSolver();
            std::vector<double> unknown_u = mySystem_u.getSolution();


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(context_.constants.k, context_.state.etaOld, context_.state.xi, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::Y, Axis::X, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_v.ThomaSolver();
            std::vector<double> unknown_v = mySystem_v.getSolution();


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(context_.constants.k, context_.state.etaOld, context_.state.xi, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::Z, Axis::X, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_w.ThomaSolver();
            std::vector<double> unknown_w = mySystem_w.getSolution();
            
            for (size_t i = 0; i < sysDimension; i++)
            {
                context_.state.eta(Axis::X, i, j, k) = unknown_u[i];
                context_.state.eta(Axis::Y, i, j, k) = unknown_v[i];
                context_.state.eta(Axis::Z, i, j, k) = unknown_w[i];
            }

        }
    }
    }
   



    // ------------------------------------------
    // Solve Zeta --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::Y;
    size_t nSystem = context_.gridPtr->Nx * context_.gridPtr->Nz; // number of linear systems to solve
    size_t sysDimension = context_.gridPtr->Ny; // dimension of linear system to solve
    // when solving Zeta we fill linsys with dyy derivatives
    // Zeta.v is then solved exploiting normal Neumann boundary conditions
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
    for (size_t i = 0; i < context_.gridPtr->Nx; i++)
    {
        for (size_t k = 0; k < context_.gridPtr->Nz; k++)
        {
            iStart = i;
            kStart = k;

            for (size_t j = 1; j < sysDimension-1; j++)
            {
                porosity = context_.constants.k(i,j,k);
                beta = 1 + (context_.timeSettings.dt * context_.constants.nu * 0.5 / porosity);
                gamma = context_.timeSettings.dt * context_.constants.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                mul = 1.0 / (context_.gridPtr->dy * context_.gridPtr->dy);
                deriv_u = (context_.state.zetaOld(Axis::X, i, j + 1, k) + context_.state.zetaOld(Axis::X, i, j - 1, k) - 2.0 * context_.state.zetaOld(Axis::X, i, j, k))*mul;
                deriv_v = (context_.state.zetaOld(Axis::Y, i, j + 1, k) + context_.state.zetaOld(Axis::Y, i, j - 1, k) - 2.0 * context_.state.zetaOld(Axis::Y, i, j, k))*mul;
                deriv_w = (context_.state.zetaOld(Axis::Z, i, j + 1, k) + context_.state.zetaOld(Axis::Z, i, j - 1, k) - 2.0 * context_.state.zetaOld(Axis::Z, i, j, k))*mul;

                rhs_u[j] = context_.state.eta(Axis::X, i,j,k) + gamma * deriv_u;
                rhs_v[j] = context_.state.eta(Axis::Y, i,j,k) + gamma * deriv_v;
                rhs_w[j] = context_.state.eta(Axis::Z, i,j,k) + gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(context_.constants.k, context_.state.zetaOld, context_.state.eta, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::X, Axis::Y, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_u.ThomaSolver();
            std::vector<double> unknown_u = mySystem_u.getSolution();


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(context_.constants.k, context_.state.zetaOld, context_.state.eta, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::Y, Axis::Y, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_v.ThomaSolver();
            std::vector<double> unknown_v = mySystem_v.getSolution();


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(context_.constants.k, context_.state.zetaOld, context_.state.eta, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::Z, Axis::Y, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_w.ThomaSolver();
            std::vector<double> unknown_w = mySystem_w.getSolution();
            
            for (size_t j = 0; j < sysDimension; j++)
            {
                context_.state.zeta(Axis::X, i, j, k) = unknown_u[j];
                context_.state.zeta(Axis::Y, i, j, k) = unknown_v[j];
                context_.state.zeta(Axis::Z, i, j, k) = unknown_w[j];
            }

        }
    }
    }





    // ------------------------------------------
    // Solve U --------------------------------
    // ------------------------------------------
    {
    normalAxis = Axis::Z;
    size_t nSystem = context_.gridPtr->Nx * context_.gridPtr->Ny; // number of linear systems to solve
    size_t sysDimension = context_.gridPtr->Nz; // dimension of linear system to solve
    // when solving U we fill linsys with dzz derivatives
    // U.w is then solved exploiting normal Neumann boundary conditions
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
    for (size_t i = 0; i < context_.gridPtr->Nx; i++)
    {
        for (size_t j = 0; j < context_.gridPtr->Ny; j++)
        {
            iStart = i;
            jStart = j;

            for (size_t k = 1; k < sysDimension-1; k++)
            {
                porosity = context_.constants.k(i,j,k);
                beta = 1 + (context_.timeSettings.dt * context_.constants.nu * 0.5 / porosity);
                gamma = context_.timeSettings.dt * context_.constants.nu * 0.5 / beta; // get gamma coefficient for each point
                // Question, if the grid is staggered, can I consider just a value of porosity for each velocity component in the grid point??

                mul = 1.0 / (context_.gridPtr->dz * context_.gridPtr->dz);
                deriv_u = (context_.state.uOld(Axis::X, i, j, k + 1) + context_.state.uOld(Axis::X, i, j, k - 1) - 2.0 * context_.state.uOld(Axis::X, i, j, k))*mul;
                deriv_v = (context_.state.uOld(Axis::Y, i, j, k + 1) + context_.state.uOld(Axis::Y, i, j, k - 1) - 2.0 * context_.state.uOld(Axis::Y, i, j, k))*mul;
                deriv_w = (context_.state.uOld(Axis::Z, i, j, k + 1) + context_.state.uOld(Axis::Z, i, j, k - 1) - 2.0 * context_.state.uOld(Axis::Z, i, j, k))*mul;

                rhs_u[k] = context_.state.zeta(Axis::X, i,j,k) + gamma * deriv_u;
                rhs_v[k] = context_.state.zeta(Axis::Y, i,j,k) + gamma * deriv_v;
                rhs_w[k] = context_.state.zeta(Axis::Z, i,j,k) + gamma * deriv_w;
            }

            mySystem_u.setRhs(rhs_u);

            mySystem_u.fillSystemVelocity(context_.constants.k, context_.state.uOld, context_.state.zeta, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::X, Axis::Z, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_u.ThomaSolver();
            std::vector<double> unknown_u = mySystem_u.getSolution();


            mySystem_v.setRhs(rhs_v);

            mySystem_v.fillSystemVelocity(context_.constants.k, context_.state.uOld, context_.state.zeta, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::Y, Axis::Z, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_v.ThomaSolver();
            std::vector<double> unknown_v = mySystem_v.getSolution();


            mySystem_w.setRhs(rhs_w);

            mySystem_w.fillSystemVelocity(context_.constants.k, context_.state.uOld, context_.state.zeta, context_.bcSettings.uBoundNew, 
                                        context_.bcSettings.uBoundOld, Axis::Z, Axis::Z, iStart, jStart, kStart, context_.constants.nu, context_.timeSettings.dt);
            mySystem_w.ThomaSolver();
            std::vector<double> unknown_w = mySystem_w.getSolution();
            
            for (size_t k = 0; k < sysDimension; k++)
            {
                context_.state.u(Axis::X, i, j, k) = unknown_u[k];
                context_.state.u(Axis::Y, i, j, k) = unknown_v[k];
                context_.state.u(Axis::Z, i, j, k) = unknown_w[k];
            }

        }
    }
    } 



    
    

}







