#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status.
# Questo previene che lo script continui se 'cmake' fallisce.
set -e
CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DPROFILING=ON"
# --- 1. Define build directory ---
BUILD_DIR="build"

# --- 2. Delete existing build directory ---
echo "--- Cleaning up old build directory: $BUILD_DIR ---"
# L'opzione -f (force) sopprime errori se la cartella non esiste
rm -rf $BUILD_DIR

# --- 3. Create new build directory ---
echo "--- Creating new build directory: $BUILD_DIR ---"
mkdir $BUILD_DIR

# --- 4. Enter the build directory ---
# Tutti i comandi successivi verranno eseguiti da dentro 'build/'
cd $BUILD_DIR

# --- 5. Run CMake (Configure/Prepare) ---
# Il '..' dice a cmake di cercare CMakeLists.txt nella cartella genitore
echo "--- Running CMake... ---"
cmake .. ${CMAKE_FLAGS}

# --- 6. Run Make (Compile/Build) ---
echo "--- Running Make... ---"
make

echo "--- Build complete! ---"
echo "You can now run your executable from within the '$BUILD_DIR' directory."
