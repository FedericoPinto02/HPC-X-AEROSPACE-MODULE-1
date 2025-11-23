#pragma once

#include <vector>
#include <memory>
#include "core/TridiagMat.hpp"
#include "numerics/LinearSys.hpp"

class SchurSequentialSolver {

public:
    SchurSequentialSolver(int globalSize, int numDomains, BoundaryType type);

    /**
     * @brief Phase 0: Pre-processing.
     * Builds the submatrices A_ii, extracts interface coefficients,
     * and assembles the Schur matrix 'S'.
     * @param A_global The GLOBAL tridiagonal matrix (N x N).
     */
    void PreProcess(const TridiagMat& A_global);

    /**
     * @brief Phase 1: Solve.
     * Solves the system A*x = f using the Schur algorithm.
     * 'PreProcess' must have been called before.
     * @param f_global The global RHS vector (length N).
     * @return The global solution vector x (length N).
     */
    std::vector<double> solve(const std::vector<double>& f_global);


private:
    // --- Decomposition variables ---
    int nGlobal;           // Total dimension N
    int num_domains;       // Number of domains P
    int num_interfaces;    // Number of interfaces (P-1)
    BoundaryType bType;

    std::vector<int> domain_sizes;  
    std::vector<int> domain_starts; 
    std::vector<int> interface_indices; 

    // --- Local and Schur systems ---
    std::vector<LinearSys> localSolvers; // P solvers (for A_ii)

    // schurSolver is created ONLY if P >= 3 (num_interfaces >= 2)
    std::unique_ptr<LinearSys> schurSolver; // <-- MODIFIED

    // schurScalarS is used ONLY if P = 2 (num_interfaces = 1)
    double schurScalarS = 0.0;

    // --- Interface matrix data ---
    // A_ie: connects u_i (interior) to u_e (interface)
    std::vector<double> A_ie_left;  // [P] 'a' (subdiagonal) coefficients that connect to the left
    std::vector<double> A_ie_right; // [P] 'c' (superdiagonal) coefficients that connect to the right

    // A_ei: connects u_e (interface) to u_i (interior)
    std::vector<double> A_ei_left;  // [P-1] 'a' (subdiagonal) coefficients that connect to the left
    std::vector<double> A_ei_right; // [P-1] 'c' (superdiagonal) coefficients that connect to the right

    // --- Private helper functions ---
    std::vector<double> getGlobalSubVector(const std::vector<double>& globalVec, 
                                           int start, int size) const;
    void setGlobalSubVector(std::vector<double>& globalVec, 
                            const std::vector<double>& subVec, 
                            int start) const;
    std::vector<double> getInterfaceVector(const std::vector<double>& globalVec) const;
    void setInterfaceVector(std::vector<double>& globalVec, 
                            const std::vector<double>& interfaceVec) const;
    std::vector<double> createSparseVector(int size, int index, double value) const;
};