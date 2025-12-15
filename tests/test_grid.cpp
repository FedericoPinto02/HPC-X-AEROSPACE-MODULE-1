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


// =========================================================================
// TEST SUITE 3: Halo & Size Verification
// =========================================================================

TEST(GridHaloTest, CalculatesSizeWithHaloCorrectly) {
    // 1. Create a Grid with known inner dimensions
    // 10x10x10 inner cells, halo=2
    size_t n_halo = 2;
    Grid grid(10, 10, 10, 1.0, 1.0, 1.0); // Uses sequential constructor

    // Manually set n_halo since sequential constructor defaults to 0 (in your cpp)
    // OR we can rely on parallel constructor if available.
    // Let's modify the struct manually to test the logic in isolation
    grid.n_halo = n_halo;

    // Inner size: 10 * 10 * 10 = 1000
    EXPECT_EQ(grid.size(), 1000);

    // Size with Halo: (10+4) * (10+4) * (10+4) = 14 * 14 * 14 = 2744
    size_t expected_dim = 10 + (2 * n_halo);
    size_t expected_total = expected_dim * expected_dim * expected_dim;

    EXPECT_EQ(grid.sizeWithHalo(), expected_total);
}

TEST(GridHaloTest, ParallelConstructorSetsHaloDefault) {
    // Use Parallel Constructor (default halo = 1)
    int argc = 0;
    char ** argv = nullptr;
    Grid grid(10, 10, 10, 1.0, 1.0, 1.0, MpiEnv(argc, argv));

    EXPECT_EQ(grid.n_halo, 1);

    // Check math roughly (if 1 proc, Nx=10. If >1 proc, Nx < 10)
    size_t expected_Nx = grid.Nx + 2;
    size_t expected_Ny = grid.Ny + 2;
    size_t expected_Nz = grid.Nz + 2;

    EXPECT_EQ(grid.sizeWithHalo(), expected_Nx * expected_Ny * expected_Nz);
}

TEST(GridHaloTest, ParallelConstructorSetsCustomHalo) {
    // Use Parallel Constructor with halo = 1
    size_t custom_halo = 1;
    int argc = 0;
    char ** argv = nullptr;
    Grid grid(10, 10, 10, 1.0, 1.0, 1.0, MpiEnv(argc, argv), custom_halo);

    EXPECT_EQ(grid.n_halo, 1);

    size_t expected_Nx = grid.Nx + 2; // + 2*1
    EXPECT_EQ(grid.sizeWithHalo(), expected_Nx * (grid.Ny + 2) * (grid.Nz + 2));
}