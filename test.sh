#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="build"
CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON"
NUM_PROCS=4      # Set the number of MPI processes for testing

# Forces OpenMPI to use a short path for temporary files
export TMPDIR=/tmp

mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR" > /dev/null

# Configure and build
cmake .. ${CMAKE_FLAGS}
make -j"$(nproc)"

popd > /dev/null

# Run all tests
echo -e "\n=== Running all tests with ${NUM_PROCS} MPI processes ==="
# - use --oversubscribe in case your machine has fewer cores than NUM_PROCS
mpirun --oversubscribe -n "${NUM_PROCS}" "${BUILD_DIR}/alltests"
