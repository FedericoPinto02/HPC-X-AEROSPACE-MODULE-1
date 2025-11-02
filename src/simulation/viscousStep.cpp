#include <algorithm>
#include <cmath>
#include <simulation/viscousStep.hpp>
#include <core/Fields.hpp>
#include <vector>
#include <numerics/derivatives.hpp>



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
}

void ViscousStep::computeG()
{
    // Ingredients list
    auto& eta = context_.state.eta;
    auto& zeta = context_.state.zeta;
    auto& u = context_.state.u;
    auto& p = context_.state.p;
    auto& k_data = context_.constants.k.getData();
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
        auto& u_data = context_.state.u(axis).getData();
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

}

void ViscousStep::closeViscousStep()
{

}


