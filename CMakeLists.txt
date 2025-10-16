cmake_minimum_required(VERSION 3.10)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)

# Include headers
include_directories(include)

# GoogleTest submodule
add_subdirectory(extern/googletest)

enable_testing()

# lib
add_library(TridiagMatrix src/core/TridiagMat.cpp)

target_include_directories(TridiagMatrix PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# executable
add_executable(testTridiag tests/test_tridiag.cpp)
target_link_libraries(testTridiag PRIVATE TridiagMatrix gtest_main)

# Add tests
include(GoogleTest)
gtest_discover_tests(testTridiag)

# === Derivatives library ===
add_library(Derivatives
    src/numerics/derivatives.cpp
    src/core/Fields.cpp   # <-- aggiungi anche questa se contiene Field::setup
)

target_include_directories(Derivatives PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# === Test for Derivatives ===
add_executable(test_derivatives tests/test_derivatives.cpp)

target_link_libraries(test_derivatives
    PRIVATE
        Derivatives
        gtest_main
)

# Register the test
gtest_discover_tests(test_derivatives)
