#pragma once
#include <vector>
#include <core/TridiagMat.hpp>
#include <core/Fields.hpp>
#include <stdexcept>
#include "simulation/SimulationContext.hpp"

enum class BoundaryType {Normal, Tangent};
class LinearSysTestFixture;

class LinearSys {

    friend class LinearSysTestFixture;
    friend class SchurSolverTestFixture;
    friend class SchurSequentialSolver;

private: 
    TridiagMat matA;
    std::vector<double> rhsC;
    std::vector<double> unknownX;
    BoundaryType boundaryType;

public:
     
    /**
     * @brief Constructor.
     * @param n Size of the linear system (n x n)
     */
    LinearSys(int n);

    /**
     * @brief Sets the right-hand side (rhsC) vector for the solver.
     * @param newRhs The vector to be copied into the internal rhsC.
     * @throws std::runtime_error if newRhs size does not match system size.
     */
    void setRhs(const std::vector<double>& newRhs);

    /**
     * @brief Gets the solution vector (unknownX) computed by the solver.
     * @return A const reference to the solution vector.
     */
    const std::vector<double>& getSolution() const;


    /**
     * @brief Fill the linear system for Pressure variables
     */
    void fillSystemPressure(const GridPtr& grid, const Axis direction);


    /**
     * @brief Fill the linear system for Velocity variables
     */
    void fillSystemVelocity(
        const SimulationData& simData, const VectorField& eta, const VectorField& xi,
        const Axis fieldComponent, const Axis derivativeDirection,
        const size_t iStart, const size_t jStart, const size_t kStart);

    /**
     * @brief Getter to access the internal matrix (needed for Schur).
     */
    const TridiagMat& getMatrix() const { return matA; }

    /**
     * @brief Getter to access the internal RHS (modified by fillSystem).
     */
    const std::vector<double>& getRhs() const { return rhsC; }

    /**
     * @brief Solve the linear system
     */
    void ThomaSolver();

};