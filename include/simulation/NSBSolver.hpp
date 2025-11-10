#ifndef NBSOLVER_NSBSOLVER_HPP
#define NBSOLVER_NSBSOLVER_HPP

#include <memory>

#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"
#include "simulation/SimulationContext.hpp"

/**
 * @brief Class responsible for solving the Navier-Stokes-Brinkman equations.
 */
class NSBSolver {
public:
    /**
     * @brief Default constructor.
     */
    NSBSolver() = default;


    /**
     * @brief Solve the Navier-Stokes-Brinkman equations.
     * @param SimulationData structure.
     */
    void solve();

private:
    SimulationData data;
    OutputSettings outputSettings;
    LoggingSettings loggingSettings;
};

#endif //NBSOLVER_NSBSOLVER_HPP
