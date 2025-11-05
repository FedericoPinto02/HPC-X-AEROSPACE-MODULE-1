#pragma once
#include <memory>
#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"
#include "simulation/SimulationData.hpp"

/**
 * @brief Class responsible for solving the Navier-Stokes-Brinkman equations.
 */
class NSBsolver {
public:
    /**
     * @brief Defalut constructor.
     */
    NSBsolver() = default;


    /**
     * @brief Solve the Navier-Stokes-Brinkman equations.
     * @param SimualtionData structure.
     */
    void solve(SimulationData& simData);

private:
   

};

