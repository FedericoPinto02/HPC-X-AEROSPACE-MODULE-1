cmake_minimum_required(VERSION 3.10)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define project-related variables
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_EXTERN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern) # for nlohmann/json

include(FetchContent)

# Include directories
include_directories(include)
include_directories(extern)  # per nlohmann/json


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

# lib
add_library(allib src/core/TridiagMat.cpp
        src/core/Fields.cpp
        src/numerics/LinearSys.cpp
        src/numerics/derivatives.cpp
        src/io/inputReader.cpp
        src/io/VTKWriter.cpp
        src/simulation/initializer.cpp
        src/simulation/viscousStep.cpp
        src/simulation/pressureStep.cpp
)
target_include_directories(allib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(allib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)


# ========================================================
# ---- TESTS ---------------------------------------------
# ========================================================
set(PROJECT_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# Option to enable/disable building tests
option(BUILD_TESTS "Build the tests" ON)
message(STATUS "BUILD_TESTS option: ${BUILD_TESTS}")

# Fetch GoogleTest source code automatically
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)

if(BUILD_TESTS)
    message(STATUS "Building tests")
    # Enable testing functionality in CMake
    enable_testing()

    # Create executable target for tests
    add_executable(alltests
            ${PROJECT_TESTS_DIR}/test_tridiag.cpp
            ${PROJECT_TESTS_DIR}/test_fields.cpp
            ${PROJECT_TESTS_DIR}/test_derivatives.cpp
            ${PROJECT_TESTS_DIR}/test_linearSys.cpp
            ${PROJECT_TESTS_DIR}/test_viscousStep.cpp
            ${PROJECT_TESTS_DIR}/test_pressureStep.cpp
    )

    # Link test executable with the libraries to be tested and GoogleTest tests runner
    target_link_libraries(alltests PRIVATE
            allib
            gtest_main
    )

    # Include GoogleTest CMake module to enable test discovery
    include(GoogleTest)
    # Automatically discover and register all GoogleTest tests in the executable
    gtest_discover_tests(alltests)

    message(STATUS "Building tests - done")
endif()