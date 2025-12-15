#include <iostream>
#include <mpi.h>
#include <vector>

#include "core/MpiEnv.hpp"
#include "numerics/ThomasSolver.hpp"

/**
 * @brief Class implementing the Schur Complement Method for solving distributed tridiagonal systems.
 *
 * The local system on each MPI process includes two interface points (left and right),
 * leading to a local system of size N = n_inner + 2, where n_inner is the number of inner points.
 * The Schur complement reduces the global system to a smaller one involving only the interface points.
 */
class SchurSolver {
private:
    static constexpr size_t NUM_LOCAL_INTERFACES = 2;
    static constexpr size_t NUM_LOCAL_SCHUR_ELEMS = 4;

    // --- MPI context ---
    int lineRank_, lineNProcs_;
    MPI_Comm lineComm_;

    // --- Local matrix coefficients (same size N) ---
    std::vector<double> a_, b_, c_;

    // --- Local chunk dimensions ---
    size_t N;        // Total local size (including interfaces at 0 and N-1)
    size_t n_inner;  // Size of the inner system (N-2)

    // --- Local Schur Complement matrix (2x2)
    // S_local = [ s00  s01 ]
    //           [ s10  s11 ]
    double s00, s01, s10, s11;

    // --- Helper solvers ---
    ThomasSolver thomas_;

public:
    /**
     * @brief Constructs the SchurSolver for a local tridiagonal system (including left and right interfaces).
     * @param env the MPI environment containing topology and local process information
     * @param a the sub-diagonal coefficients of the local tridiagonal matrix, including left and right interfaces:
     * @param b the diagonal coefficients of the local tridiagonal matrix (including left and right interfaces)
     * @param c the super-diagonal coefficients of the local tridiagonal matrix (including left and right interfaces)
     */
    SchurSolver(const MpiEnv &env,
                const Axis axis,
                const std::vector<double> &a,
                const std::vector<double> &b,
                const std::vector<double> &c)
            : a_(a), b_(b), c_(c),
              s00(0.0), s01(0.0), s10(0.0), s11(0.0), // -- temporary
              thomas_(b.size()) {
        lineComm_ = env.lineComm(axis);
        lineNProcs_ = env.lineSize(axis);
        lineRank_ = env.lineRank(axis);

        N = b.size();
        n_inner = N - 2;

        if (N < 3) { throw std::runtime_error("Local grid too small for Schur!"); }
    }

    /**
     * @brief Preprocesses the local chunk matrix (including coefficients for left and right interfaces)
     * to compute the Schur complement - i.e., the 2x2 matrix that relates the two interface points.
     * Assuming the matrix not to be time-dependent, the preprocessing can be performed only once.
     */
    void preprocess();

    /**
     * @brief Solves the full local system <code>A * u = f</code> using the Schur Complement Method.
     * @param f the local complete RHS vector (including left and right interface points)
     * @param u the local complete solution vector (to be filled)
     */
    void solve(const std::vector<double> &f, std::vector<double> &u) {
        if (u.size() != N) { u.resize(N); }
        // Phase 1: Solve for the interfaces (global communication required)
        solveInterface(f, u);
        // Phase 2: Solve for the interior (purely local)
        solveInterior(f, u);
    }

private:

    /**
     * @brief Solves the inner tridiagonal system <code>A_inner * x = rhs</code using the Thomas algorithm.
     * @param rhs the RHS vector for the inner system (size n_inner)
     * @param x the solution vector to be filled (size n_inner)
     */
    void solveInnerSystem(const std::vector<double> &rhs, std::vector<double> &x);


    /**
     * @brief Transforms the local complete RHS into a condensed one (relating to left and right interfaces)
     * - allowing to calculate the local chunk contributions for the shared interface unknowns.
     * @param f the local complete RHS (also at interface points)
     * @return the condensed RHS for the local Schur system
     */
    std::pair<double, double> condenseRHS(const std::vector<double> &f);


    /**
     * @brief Solves the interface unknowns using the Schur complement.
     * @param f the local complete RHS (including interface points)
     * @param u the local complete solution vector (to be filled at interface points)
     */
    void solveInterface(const std::vector<double> &f, std::vector<double> &u);

    /**
     * @brief Solves the 2Px2P global Schur system.
     * @param S_all the global Schur complement matrix (tridiagonal, flattened)
     *  the order is <code>[s00_proc0, s01_proc0, s10_proc0, s11_proc0, s00_proc1, ...]</code>
     * @param f_all the global condensed RHS vector
     *  the order is <code>[f_0_proc0, f_N_proc0, f_0_proc1, f_N_proc1, ...]</code>
     * @param glue_all the global coupling coefficients between processors
     *  the order is <code>[a_0_proc0, c_N-1_proc0, a_0_proc1, c_N-1_proc1, ...]</code>
     * @param u_all the global interface solution vector (to be filled);
     *  the order is <code>[u_0_proc0, u_N_proc0, u_0_proc1, u_N_proc1, ...]</code>
     */
    void solveGlobalInterfaceSystem(const std::vector<double> &S_all,
                                    const std::vector<double> &f_all,
                                    const std::vector<double> &glue_all,
                                    std::vector<double> &u_s_all);


    /**
     * @brief Solves the interior unknowns given the interface values.
     * @param f the local complete RHS vector (including interface points)
     * @param u the local complete solution vector (to be filled at interior points)
     */
    void solveInterior(const std::vector<double> &f, std::vector<double> &u);
};