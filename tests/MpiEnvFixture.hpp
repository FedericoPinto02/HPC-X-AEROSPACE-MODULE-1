#pragma once

#include <gtest/gtest.h>

#include "core/MpiEnv.hpp"

// A global accessor that all test files can use
extern MpiEnv *g_mpi;

/// GTest Environment class that manages MPI Init/Finalize automatically.
class MpiTestEnvironment : public ::testing::Environment {
private:
    int argc_;
    char **argv_;

public:
    MpiTestEnvironment(int argc, char **argv) : argc_(argc), argv_(argv) {}

    // Called once before ALL tests start
    void SetUp() override {
        // We create the MpiEnv on the heap.
        // Its constructor calls MPI_Init.
        g_mpi = new MpiEnv(argc_, argv_);
        g_mpi->setupTopology();
    }

    // Called once after ALL tests finish
    void TearDown() override {
        // Deleting it calls the destructor, which calls MPI_Finalize.
        delete g_mpi;
        g_mpi = nullptr;
    }
};
