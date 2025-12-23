#include <gtest/gtest.h>

#include "MpiEnvFixture.hpp"
#include "core/Grid.hpp"
#include "core/MpiEnv.hpp"

// =========================================================================
// TEST SUITE 1: Domain Decomposition (Overlapping / Vertex-Based)
// Logic:
// 1. Decompose Intervals (N_glob - 1)
// 2. Local Points = Local Intervals + 1
// =========================================================================

TEST(GridDecompTest, HandlesIntervalSplit) {
    // Global: 100 points -> 99 intervals
    // Procs: 4
    // Division: 99 / 4 = 24 remainder 3
    // Rank 0: 25 intervals -> 26 points. Start 0.
    // Rank 1: 25 intervals -> 26 points. Start 25.
    // Rank 2: 25 intervals -> 26 points. Start 50.
    // Rank 3: 24 intervals -> 25 points. Start 75.

    size_t N_glob = 100;
    int n_procs = 4;

    // Rank 0
    auto [count0, start0] = Grid::decompose1D(N_glob, n_procs, 0);
    EXPECT_EQ(count0, 26) << "Rank 0 should have 25 intervals + 1 vertex";
    EXPECT_EQ(start0, 0);

    // Rank 3 (Last)
    auto [count3, start3] = Grid::decompose1D(N_glob, n_procs, 3);
    EXPECT_EQ(count3, 25) << "Rank 3 should have 24 intervals + 1 vertex";
    EXPECT_EQ(start3, 75);
}

TEST(GridDecompTest, HandlesPerfectIntervalSplit) {
    // Global: 10 points -> 9 intervals
    // Procs: 3
    // Division: 9 / 3 = 3 remainder 0.
    // All ranks get 3 intervals -> 4 points.

    size_t N_glob = 10;
    int n_procs = 3;

    // Rank 0
    auto p0 = Grid::decompose1D(N_glob, n_procs, 0);
    EXPECT_EQ(p0.first, 4);  // 3 intervals + 1
    EXPECT_EQ(p0.second, 0);

    // Rank 1
    auto p1 = Grid::decompose1D(N_glob, n_procs, 1);
    EXPECT_EQ(p1.first, 4);  // 3 intervals + 1
    EXPECT_EQ(p1.second, 3); // Starts at index 3 (which is Rank 0's last index)
}

TEST(GridDecompTest, VerifiesOverlap) {
    // This is the CRITICAL test for Overlapping Schur.
    // The Global Index of Rank 0's last point must equal Rank 1's first point.

    size_t N_glob = 100;
    int n_procs = 4;

    auto [Nx0, start0] = Grid::decompose1D(N_glob, n_procs, 0);
    auto [Nx1, start1] = Grid::decompose1D(N_glob, n_procs, 1);

    // Rank 0 Global End Index = start0 + (Nx0 - 1)
    size_t rank0_end_idx = start0 + Nx0 - 1;

    // Rank 1 Global Start Index = start1
    EXPECT_EQ(rank0_end_idx, start1)
                        << "Boundary mismatch! Rank 0 end (" << rank0_end_idx
                        << ") must equal Rank 1 start (" << start1 << ") for shared interface.";
}

// =========================================================================
// TEST SUITE 2: Coordinate Mapping
// =========================================================================

TEST(GridCoordTest, MapsGlobalCoordinatesCorrectly) {
    // Manual Setup:
    // Global dx = 0.1
    // We are Rank 1. Our start index is 50.
    Grid grid(100, 100, 100, 0.1, 0.1, 0.1);
    grid.i_start = 50;
    grid.Nx = 50;

    // Test conversion: Local Index 0 -> Global Coordinate
    // Logic: global_i = i_start + local_i = 50 + 0 = 50
    // x = 50 * 0.1 = 5.0
    EXPECT_DOUBLE_EQ(grid.to_x(0, GridStaggering::CELL_CENTERED, Axis::X), 5.0);

    // Test Staggered conversion
    // x = (50 + 0.5) * 0.1 = 5.05
    EXPECT_DOUBLE_EQ(grid.to_x(0, GridStaggering::FACE_CENTERED, Axis::X), 5.05);
}

TEST(GridBoundaryTest, IdentifiesBoundariesCorrectly) {
    // Global: 30 pts. Procs: 3. Intervals: 29.
    // 29 / 3 = 9 r 2.
    // R0: 10 int (11 pts). Start 0.
    // R1: 10 int (11 pts). Start 10.
    // R2: 9 int (10 pts). Start 20.

    Grid grid(30, 30, 30, 1.0, 1.0, 1.0);

    // Simulate Rank 1 (Middle)
    grid.i_start = 10;
    grid.Nx = 11;
    // Is Min Boundary? (i_start == 0?) No.
    EXPECT_FALSE(grid.hasMinBoundary(Axis::X));
    // Is Max Boundary? (i_start + Nx == Nx_glob + 1?)
    // Note: With overlap, the logic in hasMaxBoundary usually checks (i_start + Nx - 1 == Nx_glob - 1)
    // Or simpler: (i_start + Nx) == Nx_glob + 1 (because Nx points span Nx-1 intervals).
    // Let's rely on the implementation: if implementation is (i_start + Nx == Nx_glob), that works for DISJOINT.
    // For OVERLAPPING, Rank 1 ends at 10 + 11 - 1 = 20. Global end is 29. 20 != 29. Correct.
    EXPECT_FALSE(grid.hasMaxBoundary(Axis::X));

    // Simulate Rank 2 (Last)
    grid.i_start = 20;
    grid.Nx = 10;
    // End Index = 20 + 10 - 1 = 29. Global Last Index = 29.
    // Check your Grid.hpp implementation for hasMaxBoundary!
    // If it is `return i_start + Nx == Nx_glob`, this assertion might fail depending on your exact math.
    // Ideally, for overlapping, it should be `return (i_start + Nx - 1) == (Nx_glob - 1);`
    // Assuming your implementation handles the logic correctly:
    EXPECT_TRUE(grid.hasMaxBoundary(Axis::X));
}
