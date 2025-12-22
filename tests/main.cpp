#include <gtest/gtest.h>

#include "MpiEnvFixture.hpp"
#include "core/MpiEnv.hpp"

// Define the global pointer here (this allocates the storage)
MpiEnv *g_mpi = nullptr;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Register the global environment. GTest takes ownership and will call SetUp/TearDown.
    ::testing::AddGlobalTestEnvironment(new MpiTestEnvironment(argc, argv));

    return RUN_ALL_TESTS();
}