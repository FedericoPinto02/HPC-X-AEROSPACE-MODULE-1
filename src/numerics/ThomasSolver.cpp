#include "numerics/ThomasSolver.hpp"

void ThomasSolver::solve(const std::vector<double> &a,
                         const std::vector<double> &b,
                         const std::vector<double> &c,
                         std::vector<double> &d) {

    const size_t n = d.size();
    if (n == 0) return;

    // Ensure scratch space is sufficient
    if (c_prime_.size() < n) {
        c_prime_.resize(n);
    }

    // Get raw pointers for maximum speed (avoids bounds checking)
    const double *__restrict__ pa = a.data();
    const double *__restrict__ pb = b.data();
    const double *__restrict__ pc = c.data();
    double *__restrict__ pd = d.data();
    double *__restrict__ pcp = c_prime_.data();

    // --- Forward Elimination ---
    // Step 0 (Boundary)
    double pivot = pb[0];
    if (std::abs(pivot) < EPSILON) throw std::runtime_error("ThomasSolver: Zero pivot at index 0");

    // Pre-calculate inverse to multiply instead of divide in the loop (faster)
    double inv_pivot = 1.0 / pivot;
    pcp[0] = pc[0] * inv_pivot;
    pd[0] = pd[0] * inv_pivot;

    // Step 1 to N-1
    for (size_t i = 1; i < n; ++i) {
        // m is the multiplier 1 / (b[i] - a[i] * c'[i-1])
        double m = 1.0 / (pb[i] - pa[i] * pcp[i - 1]);
        pcp[i] = pc[i] * m;
        pd[i] = (pd[i] - pa[i] * pd[i - 1]) * m;
    }

    // --- Backward Substitution ---
    // Step N-2 down to 0
    // We use a signed integer for the loop to handle the >= 0 condition safely
    for (long long i = static_cast<long long>(n) - 2; i >= 0; --i) {
        pd[i] -= pcp[i] * pd[i + 1];
    }
}