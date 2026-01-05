#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>

/**
 * @class ThomasSolver
 * @brief An efficient, specialized solver for tridiagonal linear systems.
 */
class ThomasSolver {
private:
    /// Small epsilon to avoid division by zero.
    static constexpr double EPSILON = 1e-15;

    /// Scratch buffer for modified <code>c</code> coefficients - kept as a member to avoid reallocation on every call.
    std::vector<double> c_prime_;

public:
    /**
     * @brief Constructs the solver with optional pre-allocation of the workspace.
     * @param size The expected size of the system (N). Pre-allocating prevents 
     * initial reallocations during the first solve.
     */
    explicit ThomasSolver(size_t size = 0) {
        if (size > 0) { c_prime_.resize(size); }
    }

    /**
     * @brief Solves a tridiagonal system Ax = d in-place efficiently.
     * @param a the lower diagonal coefficients (size N, a[0] is ignored)
     * @param b the main diagonal coefficients (size N)
     * @param c the upper diagonal coefficients (size N, c[N-1] is ignored)
     * @param d the RHS vector (size N), overwritten with the solution on output
     */
    void solve(const std::vector<double> &a,
               const std::vector<double> &b,
               const std::vector<double> &c,
               std::vector<double> &d);
};