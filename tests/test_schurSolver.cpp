#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include <numeric>

// Include the fixture
#include "MpiEnvFixture.hpp"
#include "numerics/TridiagMat.hpp"
#include "numerics/SchurSolver.hpp"

// --- Analytical Solution for u'' = -1, u(0)=0, u(1)=0 ---
double exact_solution(double x) {
    return 0.5 * (x - x * x);
}

TEST(SchurSolverTest, Solves1DLaplacianOverlapping) {
    // 1. Access the global MPI Environment from the fixture
    //    g_mpi is initialized in MpiTestEnvironment::SetUp
    ASSERT_NE(g_mpi, nullptr) << "MPI Environment not initialized!";
    MpiEnv& env = *g_mpi;

    int rank = env.rank();
    int size = env.size();

    // 2. Force a 1D Topology for this test
    //    We explicitly request a 1D decomposition {size, 1, 1} to ensure
    //    ranks are ordered linearly 0..N-1 along X for simple indexing logic.
    std::array<int, 3> dims = {size, 1, 1};
    std::array<int, 3> periods = {0, 0, 0};
    try {
        env.setupTopology(dims, periods);
    } catch (const std::exception& e) {
        // If setup fails (e.g. topology locked), we might need to skip or fail.
        // For this test logic to work mathematically, we really need 1D.
        GTEST_SKIP() << "Could not enforce 1D topology: " << e.what();
    }

    // 3. Grid Parameters (Overlapping)
    //    In an overlapping decomposition, neighboring ranks share 1 point.
    //    N_global = Sum(N_local - 1) + 1
    int points_per_proc = 20;
    int local_N = points_per_proc;

    // Calculate global size based on overlap
    // rank 0 has [0..19], rank 1 has [19..38], etc.
    // Total intervals = size * (local_N - 1)
    // Total points = Total intervals + 1
    int global_N = size * (local_N - 1) + 1;

    double L = 1.0;
    double dx = L / (static_cast<double>(global_N) - 1.0);

    // 4. Determine Global Start Index
    //    Rank i starts at index i * (local_N - 1)
    int global_start_idx = rank * (local_N - 1);

    // 5. Build Matrix using TridiagMat
    TridiagMat matrix(local_N);

    // Get references to diagonals
    std::vector<double>& sub = matrix.getDiag(-1);
    std::vector<double>& main = matrix.getDiag(0);
    std::vector<double>& sup = matrix.getDiag(1);

    // --- FILL MATRIX (Laplacian Stencil) ---
    // Stencil: [-1, 2, -1] / dx^2
    // We put the 1/dx^2 scaling into the RHS for simplicity, or keep coeffs raw.
    // Here: Matrix * u = dx^2 * f(x).
    std::fill(sub.begin(), sub.end(), -1.0);
    std::fill(main.begin(), main.end(), 2.0);
    std::fill(sup.begin(), sup.end(), -1.0);

    // --- ISOLATE LOCAL CHUNK ---
    // The solver must treat the local matrix as a standalone block.
    // Connections to ghosts are handled by the Schur complement Schur matrix S.
    sub[0] = 0.0;            // No dependency on left ghost
    sup[local_N - 1] = 0.0;  // No dependency on right ghost

    // 6. Build RHS
    //    -u'' = 1  =>  RHS = 1 * dx^2
    std::vector<double> rhs(local_N, dx * dx);

    // 7. Apply Boundary Conditions & Overlap Scaling
    //    Rule: For shared nodes (internal interfaces), we scale the operator
    //    (diagonal and RHS) by 0.5 so the sum across ranks equals the full operator.

    // --- LEFT END of Local Chunk ---
    if (rank == 0) {
        // Physical Boundary: Dirichlet u=0
        main[0] = 1.0;
        sup[0] = 0.0; // Cut connection to right neighbor (row becomes identity)
        rhs[0] = 0.0;
    } else {
        // Internal Interface (Shared with Rank i-1)
        // Scale contributions
        main[0] *= 0.5;
        rhs[0] *= 0.5;
    }

    // --- RIGHT END of Local Chunk ---
    if (rank == size - 1) {
        // Physical Boundary: Dirichlet u=0
        main[local_N - 1] = 1.0;
        sub[local_N - 1] = 0.0; // Cut connection to left neighbor
        rhs[local_N - 1] = 0.0;
    } else {
        // Internal Interface (Shared with Rank i+1)
        // Scale contributions
        main[local_N - 1] *= 0.5;
        rhs[local_N - 1] *= 0.5;
    }

    // 8. Instantiate and Run Schur Solver
    //    We use Axis::X because we set up a 1D topology along X.
    SchurSolver solver(env, Axis::X, matrix);

    solver.preprocess(); // Compute S matrix and local Green's functions

    std::vector<double> u(local_N);
    solver.solve(rhs, u);

    // 9. Verification
    double max_local_error = 0.0;
    for (int i = 0; i < local_N; ++i) {
        // Map local index 'i' to global physical coordinate 'x'
        int global_i = global_start_idx + i;
        double x = global_i * dx;
        double u_exact = exact_solution(x);

        double err = std::abs(u[i] - u_exact);
        if (err > max_local_error) max_local_error = err;
    }

    // Gather max error from all ranks
    double global_max_error = 0.0;
    MPI_Allreduce(&max_local_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0) {
        // Optional debug output
        // std::cout << "[  INFO ] Global N: " << global_N << ", dx: " << dx << std::endl;
        // std::cout << "[  INFO ] Max Error: " << global_max_error << std::endl;
    }

    // Finite Difference Accuracy Check
    // Error should be O(dx^2). For ~60-80 points, dx ~ 0.015, error ~ 1e-4.
    // We set a safe upper bound.
    EXPECT_LT(global_max_error, 1e-3);
}