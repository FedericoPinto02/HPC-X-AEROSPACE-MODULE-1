cmake_minimum_required(VERSION 3.10)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)
include_directories(extern)

# Libraries
add_library(TridiagMatrix src/core/TridiagMat.cpp)
target_include_directories(TridiagMatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Fields src/core/Fields.cpp)
target_include_directories(Fields PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Derivatives src/numerics/derivatives.cpp src/core/Fields.cpp)
target_include_directories(Derivatives PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(InputReader src/io/inputReader.cpp)
target_include_directories(InputReader PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include 
    ${CMAKE_CURRENT_SOURCE_DIR}/extern
)

add_library(Initializer src/simulation/initializer.cpp)
target_include_directories(Initializer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Initializer PUBLIC Fields InputReader)

add_library(VTKWriter src/io/VTKWriter.cpp)
target_include_directories(VTKWriter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Executable (main)
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE
    InputReader
    Initializer
    Fields
    VTKWriter
)

# GoogleTest
add_subdirectory(extern/googletest)

enable_testing()

# Libraries for tests
add_library(allib
    src/core/TridiagMat.cpp
    src/core/Fields.cpp
    src/numerics/derivatives.cpp
)
target_include_directories(allib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Executable (tests)
add_executable(alltests
    tests/test_tridiag.cpp
    tests/test_fields.cpp
    tests/test_derivatives.cpp
)

target_link_libraries(alltests PRIVATE allib gtest_main)

# Register tests
include(GoogleTest)
gtest_discover_tests(alltests)
