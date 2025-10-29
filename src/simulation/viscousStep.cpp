#include <algorithm>
#include <cmath>
#include <simulation/viscousStep.hpp>
#include <core/Fields.hpp>
#include <vector>



ViscousStep::ViscousStep(SimulationContext& ctx) : context_(ctx)
{
}

void ViscousStep::run()
{
    computeG();
    computeXi();
    closeViscousStep();
}

void ViscousStep::computeG()
{
    // Ingredients list
    auto& f = context_.constants.f; 
    auto& nu = context_.constants.nu;
    auto& k = context_.constants.k;
    auto& eta = context_.state.eta;
    auto& zeta = context_.state.zeta;
    auto& u = context_.state.u;
    auto& p = context_.state.p;

    // Recepie
    // g = f -grad(p) -nu/2/k*u +nu/2*(dxx eta + dyy zeta + dzz u)

    // Let me cook
    VectorField g;


}

void ViscousStep::computeXi()
{

}

void ViscousStep::closeViscousStep()
{

}


