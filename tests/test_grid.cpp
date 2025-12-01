#include <gtest/gtest.h>
#include "core/Grid.hpp"
#include "core/MpiEnv.hpp"

// =========================================================================
// TEST SUITE 1: Domain Decomposition (Static Math Logic)
// =========================================================================

TEST(GridDecompTest, HandlesPerfectSplit) {
    size_t N_glob = 100;
    int n_procs = 4;

    // Rank 0
    auto [count0, start0] = Grid::decompose1D(N_glob, n_procs, 0);
    EXPECT_EQ(count0, 25);
    EXPECT_EQ(start0, 0);

    // Rank 3 (Last)
    auto [count3, start3] = Grid::decompose1D(N_glob, n_procs, 3);
    EXPECT_EQ(count3, 25);
    EXPECT_EQ(start3, 75);
}

TEST(GridDecompTest, HandlesRemaindersCorrectly) {
    // Case: 10 points, 3 procs. Remainder = 1.
    // Rank 0 should take the extra point.
    size_t N_glob = 10;
    int n_procs = 3;

    // Rank 0: Base(3) + 1 = 4. Start 0.
    auto p0 = Grid::decompose1D(N_glob, n_procs, 0);
    EXPECT_EQ(p0.first, 4);
    EXPECT_EQ(p0.second, 0);

    // Rank 1: Base(3). Start 4 (0 + 4).
    auto p1 = Grid::decompose1D(N_glob, n_procs, 1);
    EXPECT_EQ(p1.first, 3);
    EXPECT_EQ(p1.second, 4);
}

// =========================================================================
// TEST SUITE 2: Coordinate Mapping
// =========================================================================

TEST(GridCoordTest, MapsGlobalCoordinatesCorrectly) {
    // Use Sequential constructor to manually set up a state
    // mimicking a "middle" block
    Grid grid(100, 100, 100, 0.1, 0.1, 0.1);

    // Manually force offsets (simulating we are Rank 1 in a 2-proc split)
    grid.i_start = 50;
    grid.Nx = 50;

    // Test conversion: Local Index 0 -> Global Coordinate
    // Formula: x = (i_start + i) * dx => (50 + 0) * 0.1 = 5.0
    EXPECT_DOUBLE_EQ(grid.to_x(0, GridStaggering::CELL_CENTERED, Axis::X), 5.0);

    // Test Staggered conversion
    // Formula: x = (i_start + i + 0.5) * dx => (50 + 0 + 0.5) * 0.1 = 5.05
    EXPECT_DOUBLE_EQ(grid.to_x(0, GridStaggering::FACE_CENTERED, Axis::X), 5.05);
}

TEST(GridBoundaryTest, IdentifiesBoundariesCorrectly) {
    // Setup: We are the middle block in a 3-block split (Offsets: 0, 10, 20)
    Grid grid(30, 30, 30, 1.0, 1.0, 1.0);
    grid.Nx = 10;
    grid.i_start = 10; // We are in the middle

    // Should NOT be min or max boundary
    EXPECT_FALSE(grid.hasMinBoundary(Axis::X));
    EXPECT_FALSE(grid.hasMaxBoundary(Axis::X));

    // Setup: We are the LAST block
    grid.i_start = 20;
    EXPECT_TRUE(grid.hasMaxBoundary(Axis::X));
}
