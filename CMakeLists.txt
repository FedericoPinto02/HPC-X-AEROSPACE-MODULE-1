cmake_minimum_required(VERSION 3.10)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
include_directories(include)
include_directories(extern)  # per nlohmann/json


# GoogleTest submodule
add_subdirectory(extern/googletest)

enable_testing()

# lib
add_library(allib src/core/TridiagMat.cpp 
                  src/core/Fields.cpp
                  src/numerics/derivatives.cpp
)

target_include_directories(allib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


# executable
add_executable(alltests tests/test_tridiag.cpp
                        tests/test_fields.cpp
                        tests/test_derivatives.cpp
)


target_link_libraries(alltests PRIVATE allib gtest_main)




# Add tests
include(GoogleTest)
gtest_discover_tests(alltests)
