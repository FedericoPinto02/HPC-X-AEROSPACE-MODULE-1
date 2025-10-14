cmake_minimum_required(VERSION 3.10)
project(NBSolver)

set(CMAKE_CXX_STANDARD 17)

# Include headers
include_directories(include)

# Add GoogleTest submodule
add_subdirectory(external/googletest)

enable_testing()


# Add tests
include(GoogleTest)
gtest_discover_tests(test_grayscale)