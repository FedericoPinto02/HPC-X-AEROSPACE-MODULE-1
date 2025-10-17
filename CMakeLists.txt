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
add_library(Fields src/core/Fields.cpp)

target_include_directories(TridiagMatrix PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_include_directories(Fields PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# executable
add_executable(testTridiag tests/test_tridiag.cpp)
target_link_libraries(testTridiag PRIVATE TridiagMatrix gtest_main)

add_executable(testFields tests/test_fields.cpp)
target_link_libraries(testFields PRIVATE Fields gtest_main)

# Add tests
include(GoogleTest)
gtest_discover_tests(testTridiag)
gtest_discover_tests(testFields)