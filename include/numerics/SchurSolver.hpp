#pragma once

#include <iostream>
#include <vector>

#include "core/MpiEnv.hpp"
#include "core/TridiagMat.hpp"
#include "numerics/ThomasSolver.hpp"

/**
 * @brief Class implementing the Schur Complement Method for solving distributed tridiagonal systems.
 *
 * The local system on each MPI process includes two interface points (left and right),
 * leading to a local system of size <code>N = n_inner + 2</code>, where n_inner is the number of inner points.
 * The Schur complement reduces the global system to a smaller one involving only the interface points.
 */
class SchurSolver {
private:
    static constexpr size_t NUM_LOCAL_INTERFACES = 2;
    static constexpr size_t NUM_LOCAL_SCHUR_ELEMS = 4;

    // --- MPI context ---
    int lineRank_, lineNProcs_;
    MPI_Comm lineComm_;

    // --- Local chunk dimensions ---
    size_t N;        // Total local size (including interfaces at 0 and N-1)
    size_t n_inner;  // Size of the inner system (N-2)
    size_t M;        // Size of shared interfaces (lineNProcs_ + 1) i.e., size of the global Schur system

    // --- Local matrix coefficients (same size N) ---
    std::vector<double> a_, b_, c_;
    // --- Global (rank 0 only) Schur matrix coefficients [(nProc+1)x(nProc+1)] ---
    std::vector<double> Sa_glob_, Sb_glob_, Sc_glob_;

    // --- Scratchpads to avoid reallocations ---
    std::vector<double> a_in_, b_in_, c_in_;    // scratchpad for A_ii
    std::vector<double> f_in_, y_;              // scratchpads for condenseRHS() and solveInterior() rhss / solutions

    // --- Helper solvers ---
    ThomasSolver thomas_in_;
    ThomasSolver thomas_s_glob_;

public:
    /**
     * @brief Constructs the SchurSolver for a local tridiagonal system (including left and right interfaces).
     * @param env the MPI environment containing topology and local process information
     * @param a the sub-diagonal coefficients of the local tridiagonal matrix, including left and right interfaces:
     * @param b the diagonal coefficients of the local tridiagonal matrix (including left and right interfaces)
     * @param c the super-diagonal coefficients of the local tridiagonal matrix (including left and right interfaces)
     */
    SchurSolver(const MpiEnv &env, const Axis axis, const TridiagMat &matrix);

    /**
     * @brief Preprocesses the local chunk matrix (including coefficients for left and right interfaces)
     * to compute the local Schur complements (i.e., the 2x2 matrix that relates the two interface points)
     * and gather them into the global Schur matrix (on the line master only).
     * If the matrix is neither time-dependent nor space-dependent, the preprocessing can be performed only once.
     */
    void preprocess();

    /**
     * @brief Solves the full local system <code>A * u = f</code> using the Schur Complement Method.
     * @param f the local complete RHS vector (including left and right interface points)
     * @param u the local complete solution vector (to be filled)
     */
    void solve(const std::vector<double> &f, std::vector<double> &u);

private:

    /**
     * @brief Solves the inner tridiagonal system <code>A_inner * x = rhs</code using the Thomas algorithm.
     * @param rhs the RHS vector for the inner system (size n_inner)
     * @param x the solution vector to be filled (size n_inner)
     */
    void inline solveInnerSystem(const std::vector<double> &rhs, std::vector<double> &x) {
        // Copy RHS into x (ThomasSolver solves in-place)
        x = rhs;
        thomas_in_.solve(a_in_, b_in_, c_in_, x);
    }


    /**
     * @brief Transforms the local complete RHS into a condensed one (relating to left and right interfaces)
     * - allowing to calculate the local chunk contributions for the shared interface unknowns.
     * @param f the local complete RHS (also at interface points)
     * @return the condensed RHS for the local Schur system
     */
    std::pair<double, double> condenseRHS(const std::vector<double> &f);


    /**
     * @brief Calculates the unknowns on the interfaces by condensing the local RHSs following the Schur method,
     *  gathering the condensed RHSs from all processes, and solving against the global Schur matrix.
     * @param f the local complete RHS (including interface points)
     * @param u the local complete solution vector (to be filled at interface points)
     */
    void solveInterface(const std::vector<double> &f, std::vector<double> &u);

    /**
     * @brief Calculates the unknowns on the interfaces by solving the ~PxP global Schur system.
     * @param f_all the global condensed RHS vector
     *  the order is <code>[f_0_proc0, f_N_proc0, f_0_proc1, f_N_proc1, ...]</code>
     * @param u_all the global interface solution vector (to be filled);
     *  the order is <code>[u_0_proc0, u_N_proc0, u_0_proc1, u_N_proc1, ...]</code>
     */
    void solveGlobalInterfaceSystem(const std::vector<double> &f_all,
                                    std::vector<double> &u_s_all);


    /**
     * @brief Solves the interior unknowns given the interface values.
     * @param f the local complete RHS vector (including interface points)
     * @param u the local complete solution vector (to be filled at interior points)
     */
    void solveInterior(const std::vector<double> &f, std::vector<double> &u);
};