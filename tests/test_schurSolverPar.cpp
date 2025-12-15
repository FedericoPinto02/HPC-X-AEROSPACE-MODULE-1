#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <array>

#include "core/MpiEnv.hpp"
#include "numerics/SchurSolver.hpp"

// --- Analytical Solution for u'' = -1, u(0)=0, u(1)=0 ---
double exact_solution(double x) {
    return 0.5 * (x - x * x);
}

TEST(SchurSolverTest, Solves1DLaplacianCorrectly) {
    // 1. Initialize MPI Environment
    // We pass dummy args since GTest main usually handles real args or we assume MPI is up.
    // Note: This relies on MpiEnv checking if MPI is already initialized.
    int argc = 0;
    char **argv = nullptr;
    MpiEnv env(argc, argv);

    int rank = env.rank();
    int size = env.size();

    // 2. Setup Topology
    // We force a 1D decomposition along the X-axis for this test.
    // dims = {size, 1, 1} means all processors are in a single line along X.
    std::array<int, 3> dims = {size, 1, 1};
    std::array<int, 3> periods = {0, 0, 0};

    // Check if enough procs are available to allow testing (SchurSolver needs N >= 3 locally)
    // We want a global grid large enough.
    env.setupTopology(dims, periods);

    // 3. Define Grid Parameters
    // We want a decent resolution. Let's say 20 points per processor.
    int points_per_proc = 20;
    int global_N = points_per_proc * size;
    double L = 1.0;
    double dx = L / (global_N - 1); // Grid spacing

    // 4. Determine Local Grid Chunk
    // Since we divided evenly (points_per_proc * size), math is simple.
    int local_N = points_per_proc;
    int global_start_idx = rank * local_N;

    // 5. Construct Local Tridiagonal System
    // Discretization of -u'' = 1  =>  (-u_L + 2u_C - u_R) / dx^2 = 1
    // Linear system form: -u_L + 2u_C - u_R = dx^2
    std::vector<double> a(local_N, -1.0); // Lower
    std::vector<double> b(local_N, 2.0); // Main
    std::vector<double> c(local_N, -1.0); // Upper
    std::vector<double> rhs(local_N, dx * dx);

    // 6. Apply Boundary Conditions (Dirichlet u=0)
    // IMPORTANT: The Schur solver expects the system for interfaces to be valid.
    // We implement Dirichlet BCs by modifying the matrix rows at the physical boundaries.

    // Left Boundary (Global Node 0)
    if (rank == 0) {
        // Equation: 1 * u_0 = 0
        b[0] = 1.0;
        c[0] = 0.0;
        // a[0] is theoretically -1 connecting to "left of 0", but boundary cuts it. Set to 0.
        a[0] = 0.0;
        rhs[0] = 0.0;
    }

    // Right Boundary (Global Node N-1)
    if (rank == size - 1) {
        // Equation: 1 * u_N = 0
        b[local_N - 1] = 1.0;
        // c[local_N - 1] connecting to right is 0.
        c[local_N - 1] = 0.0;
        a[local_N - 1] = 0.0;
        rhs[local_N - 1] = 0.0;
    }

    // 7. Instantiate and Run Schur Solver
    // We use Axis::X because we split the topology along X.
    SchurSolver solver(env, Axis::X, a, b, c);

    solver.preprocess(); // Compute S matrix

    std::vector<double> u(local_N);
    solver.solve(rhs, u);

    // 8. Verification
    double max_local_error = 0.0;
    for (int i = 0; i < local_N; ++i) {
        int global_i = global_start_idx + i;
        double x = global_i * dx;
        double u_analytical = exact_solution(x);

        double error = std::abs(u[i] - u_analytical);
        if (error > max_local_error) {
            max_local_error = error;
        }
    }

    // Reduce error across all processors to get global max error
    double global_max_error = 0.0;
    MPI_Allreduce(&max_local_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "[  INFO ] Global Max Error: " << global_max_error << std::endl;
        std::cout << "[  INFO ] Grid Spacing dx: " << dx << std::endl;
    }

    // 9. Assertion
    // Finite Difference error is O(dx^2).
    // For dx approx 0.01 (if N=100), error should be small (approx 1e-4 range).
    // We verify it's below a sensible threshold.
    EXPECT_LT(global_max_error, 1e-3);
}
