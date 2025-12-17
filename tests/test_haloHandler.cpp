#include <gtest/gtest.h>
#include <cmath>

#include "MpiEnvFixture.hpp"  // Assumes g_mpi is defined here
#include "core/Grid.hpp"
#include "core/Fields.hpp"
#include "core/HaloHandler.hpp" // or "core/HaloCommunicator.hpp" depending on your filename

class HaloHandlerTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Ensure g_mpi is initialized
        ASSERT_NE(g_mpi, nullptr);
    }
};

/**
 * @brief Test X-Direction Exchange with Overlapping Logic.
 * * SCENARIO:
 * Two processes in X (Rank 0 and Rank 1).
 * They share a single interface node.
 * * ASSERTIONS:
 * 1. Rank 0's Right Halo must receive Rank 1's FIRST INTERNAL node (index 1), NOT the boundary (index 0).
 * 2. Rank 1's Left Halo must receive Rank 0's LAST INTERNAL node (index Nx-2), NOT the boundary (index Nx-1).
 */
TEST_F(HaloHandlerTest, ExchangeX_SkipsSharedInterface) {
    if (g_mpi->size() < 2) {
        GTEST_SKIP() << "Test requires at least 2 MPI processes.";
    }

    // 1. Setup Grid: 10 Global points, 2 Procs -> Split in X
    // Rank 0: [0..5] (size 6)
    // Rank 1: [5..9] (size 5)
    // Shared Node Global Index: 5
    size_t N_glob = 10;
    Grid grid(N_glob, 3, 3, 1.0, 1.0, 1.0, *g_mpi);

    // Ensure we are split in X for this test logic to hold
    if (g_mpi->dims()[0] < 2) {
        GTEST_SKIP() << "MPI Topology is not split in X direction.";
    }

    // 2. Setup Field with Global Linear Data: f(x) = global_i
    // This makes verification easy: value == global index.
    Field field;
    field.setup(std::make_shared<Grid>(grid));

    // Populate with global indices
    auto &data = field.getData();
    for (long k = 0; k < grid.Nz; ++k) {
        for (long j = 0; j < grid.Ny; ++j) {
            for (long i = 0; i < grid.Nx; ++i) {
                // Calculate global X index
                long global_i = grid.i_start + i;
                field(i, j, k) = static_cast<double>(global_i);
            }
        }
    }

    // 3. Perform Exchange
    HaloHandler handler(*g_mpi);
    handler.exchangeX(field);

    // 4. Verification
    int rank = g_mpi->rank();
    int coords[3];
    MPI_Cart_coords(g_mpi->cartComm(), rank, 3, coords);

    // Only check if we are on the specific boundary ranks for X
    bool isLeftRank = (coords[0] == 0);
    bool isRightRank = (coords[0] == 1);

    // Check Z=1, Y=1 (middle of the small chunk)
    long j = 1;
    long k = 1;

    if (isLeftRank) {
        // I am Rank 0. My Right Boundary is Global 5.
        // My Right Halo (local Nx) should receive Global 6.
        // (Because Global 5 is shared, so we skip it and grab Global 6).

        // field(Nx, j, k) accesses the right halo
        double haloVal = field(grid.Nx, j, k);
        EXPECT_NEAR(haloVal, 6.0, 1e-9)
                            << "Rank 0 Right Halo should be Global 6 (Rank 1's internal node), but got " << haloVal;
    }

    if (isRightRank) {
        // I am Rank 1. My Left Boundary is Global 5.
        // My Left Halo (local -1) should receive Global 4.
        // (Because Global 5 is shared, so we skip it and grab Global 4).

        // field(-1, j, k) accesses the left halo
        double haloVal = field.valueWithOffset(0, j, k, Axis::X, -1);
        EXPECT_NEAR(haloVal, 4.0, 1e-9)
                            << "Rank 1 Left Halo should be Global 4 (Rank 0's internal node), but got " << haloVal;
    }
}

/**
 * @brief Test Corner Propagation (Diagonal Exchange).
 * * SCENARIO:
 * 2x2 Process Grid in X/Y.
 * We rely on the sequential nature of exchange (X then Y) to fill the corner halos.
 * * LOGIC:
 * 1. Rank (0,0) sends data to Rank (1,0) [Right Neighbor].
 * 2. Rank (1,0) receives it into its Left Halo.
 * 3. Rank (1,0) performs Y-exchange. It sends its Left Halo (which contains data from 0,0) to Rank (1,1) [Top Neighbor].
 * 4. Rank (1,1) receives it into its Bottom-Left Corner Halo.
 */
TEST_F(HaloHandlerTest, SequentialExchange_FillsCorners) {
    if (g_mpi->size() < 4) {
        GTEST_SKIP() << "Test requires at least 4 MPI processes (2x2 split).";
    }

    // 1. Setup Grid 2x2 split
    Grid grid(10, 10, 3, 1.0, 1.0, 1.0, *g_mpi);
    if (g_mpi->dims()[0] < 2 || g_mpi->dims()[1] < 2) {
        GTEST_SKIP() << "MPI Topology is not split in X and Y.";
    }

    // 2. Setup Field: f(x,y) = x + y
    Field field;
    field.setup(std::make_shared<Grid>(grid));

    // Reset to -999 to detect uninitialized halos
    field.reset(-999.0);

    for (long k = 0; k < grid.Nz; ++k) {
        for (long j = 0; j < grid.Ny; ++j) {
            for (long i = 0; i < grid.Nx; ++i) {
                double val = (grid.i_start + i) + (grid.j_start + j);
                field(i, j, k) = val;
            }
        }
    }

    // 3. Full Exchange (X -> Y -> Z)
    HaloHandler handler(*g_mpi);
    handler.exchange(field);

    // 4. Verify Corner on Rank (1,1) [Top-Right Proc]
    int rank = g_mpi->rank();
    int coords[3];
    MPI_Cart_coords(g_mpi->cartComm(), rank, 3, coords);

    if (coords[0] == 1 && coords[1] == 1) {
        // I am the Top-Right processor.
        // My Left-Bottom Corner Ghost (-1, -1) should come from Rank (0,0).
        // Rank (0,0) ends at Global X=5, Y=5. (Assuming balanced split of 10).
        // Actually, let's look at the Overlap logic:
        // Rank 00 ends at 5. Rank 11 starts at 5.
        // Shared node is (5,5).
        // My (-1, -1) ghost should be (Global 4, Global 4).

        double expected = 4.0 + 4.0; // 8.0

        // Access (-1, -1, 0)
        // We use valueWithOffset twice or manually via idx if exposed,
        // but let's assume valueWithOffset chaining or raw access logic:
        // field.valueWithOffset(0, 0, 0, Axis::X, -1) gives (-1, 0, 0)
        // To get (-1, -1, 0), we need to check internal memory layout or rely on implementation correctness.
        // Let's use raw data index:

        const size_t H = grid.n_halo;
        const size_t Nx_tot = grid.Nx + 2 * H;
        const size_t Ny_tot = grid.Ny + 2 * H;

        // Index of (-1, -1, 0) corresponds to (H-1, H-1, H) in raw buffer
        // because physical data starts at (H, H, H).
        size_t idx = ((0 + H) * Ny_tot + (H - 1)) * Nx_tot + (H - 1);

        double val = field.getData()[idx];

        EXPECT_NEAR(val, expected, 1e-9)
                            << "Corner ghost (-1, -1) on Rank (1,1) incorrect. "
                            << "Should be derived from Rank(0,0) internal node.";
    }
}

TEST_F(HaloHandlerTest, VectorField_Exchange_WorksForAllComponents) {
    if (g_mpi->size() < 2) GTEST_SKIP();

    Grid grid(10, 3, 3, 1.0, 1.0, 1.0, *g_mpi);
    VectorField vField;
    vField.setup(std::make_shared<Grid>(grid));

    // Initialize: U=1, V=2, W=3
    vField(Axis::X).reset(1.0);
    vField(Axis::Y).reset(2.0);
    vField(Axis::Z).reset(3.0);

    // Perturb boundaries to 0 to ensure we actually see halos coming in
    // (If neighbor sends 1.0, and I have 0.0, successful exchange = 1.0)
    // Actually, reset() sets halos to 1.0 too.
    // Let's set inner domain to ID, and check halos.

    // Reset everything to -1
    vField(Axis::X).reset(-1.0);
    vField(Axis::Y).reset(-2.0);
    vField(Axis::Z).reset(-3.0);

    // Set Inner to correct ID
    auto setInner = [&](Field &f, double val) {
        for (long k = 0; k < grid.Nz; ++k)
            for (long j = 0; j < grid.Ny; ++j)
                for (long i = 0; i < grid.Nx; ++i)
                    f(i, j, k) = val;
    };

    setInner(vField(Axis::X), 1.0);
    setInner(vField(Axis::Y), 2.0);
    setInner(vField(Axis::Z), 3.0);

    // Exchange
    HaloHandler handler(*g_mpi);
    handler.exchange(vField);

    // Check X-Halo of X-Component. Should be 1.0 (received from neighbor)
    // field(Nx, 0, 0)
    if (g_mpi->dims()[0] > 1) {
        double val = vField(Axis::X)(grid.Nx, 0, 0);
        // If I have a right neighbor, I expect 1.0. If boundary, it remains -1.0 (or boundary logic).
        // Since we didn't implement BCs in HaloHandler (only MPI), boundary halos stay untouched.

        // Only check if we are internal or have neighbor
        int coords[3];
        MPI_Cart_coords(g_mpi->cartComm(), g_mpi->rank(), 3, coords);
        if (coords[0] < g_mpi->dims()[0] - 1) {
            EXPECT_EQ(val, 1.0);
        }
    }
}