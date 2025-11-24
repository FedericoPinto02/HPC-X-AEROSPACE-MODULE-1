#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="build"

CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=OFF -DPROFILING=ON"
mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR" > /dev/null

# Configure and build
cmake .. ${CMAKE_FLAGS}
make -j"$(nproc)"

# Run the main program
echo -e "\n=== Running main program ==="
./main

popd > /dev/null
