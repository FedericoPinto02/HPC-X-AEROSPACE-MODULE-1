cmake_minimum_required(VERSION 3.15)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_RELEASE
        "-O3 -DNDEBUG -march=native -flto -funroll-loops -fno-math-errno")
set(CMAKE_EXE_LINKER_FLAGS "-flto")

# Find MPI
find_package(MPI REQUIRED)
message(STATUS "MPI Found: ${MPI_CXX_COMPILER}")

include(FetchContent)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ========================================================
# ---- SOLVER LIBRARY ------------------------------------
# ========================================================
# Define project-related variables
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_EXTERN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern) # for nlohmann/json

# todo - use it when all files actually compile ! (allow WIP files not to be compiled and break everything :D)
# file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/*.cpp)
# list(FILTER SOURCE_FILES EXCLUDE REGEX ".*main.cpp$") # exclude the main.cpp
# file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS ${PROJECT_INCLUDE_DIR}/*.hpp)

# Add owned src and headers in the library
add_library(NSBSolver
        # ${SOURCE_FILES} ${HEADER_FILES}   todo - use it when all files actually compile !
        # Explicit header files
        include/core/Fields.hpp
        include/core/Grid.hpp
        include/core/MpiEnv.hpp
        include/core/TridiagMat.hpp
        include/io/inputReader.hpp
        include/io/VTKWriter.hpp
        include/io/logWriter.hpp
        include/numerics/derivatives.hpp
        include/numerics/LinearSys.hpp
        include/numerics/SchurSequentialSolver.hpp
        include/numerics/SchurSolver.hpp
        include/numerics/ThomasSolver.hpp
        include/simulation/initializer.hpp
        include/simulation/pressureStep.hpp
        include/simulation/viscousStep.hpp
        include/simulation/NSBSolver.hpp
        # Explicit source files
        src/core/Fields.cpp
        src/core/Grid.cpp
        src/core/MpiEnv.cpp
        src/core/TridiagMat.cpp
        src/io/inputReader.cpp
        src/io/VTKWriter.cpp
        src/io/logWriter.cpp
        # -- src/numerics/derivatives.cpp use this source file to compile a non-optimized version
        src/numerics/derivativesOpt.cpp
        #src/numerics/derivatives.cpp
        src/numerics/LinearSys.cpp
        src/numerics/SchurSequentialSolver.cpp
        src/numerics/SchurSolver.cpp
        src/numerics/ThomasSolver.cpp
        src/simulation/initializer.cpp
        src/simulation/pressureStep.cpp
        src/simulation/viscousStep.cpp
        src/simulation/NSBSolver.cpp

)
# Include owned headers and external headers.
target_include_directories(NSBSolver PUBLIC
        ${PROJECT_INCLUDE_DIR}
        ${PROJECT_EXTERN_DIR}
)
# Link MPI to the solver library
target_link_libraries(NSBSolver PUBLIC MPI::MPI_CXX)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ========================================================
# ---- APPLICATIONS --------------------------------------
# ========================================================
add_executable(main ${PROJECT_SOURCE_DIR}/main.cpp)
target_link_libraries(main PRIVATE NSBSolver)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ========================================================
# ---- TESTS ---------------------------------------------
# ========================================================
# Define test-related variables
set(PROJECT_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# Option to enable/disable building tests
option(BUILD_TESTS "Build the tests" ON)
message(STATUS "BUILD_TESTS option: ${BUILD_TESTS}")

# Fetch GoogleTest source code automatically
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)

if (BUILD_TESTS)
    message(STATUS "Building tests")
    # Enable testing functionality in CMake
    enable_testing()

    # Create executable target for tests
    add_executable(alltests
            ${PROJECT_TESTS_DIR}/test_derivatives.cpp
            ${PROJECT_TESTS_DIR}/test_grid.cpp
            ${PROJECT_TESTS_DIR}/test_fields.cpp
            ${PROJECT_TESTS_DIR}/test_pressureStep.cpp
            ${PROJECT_TESTS_DIR}/test_viscousStep.cpp
            ${PROJECT_TESTS_DIR}/test_schurSolver.cpp
            ${PROJECT_TESTS_DIR}/test_tridiag.cpp
            #${PROJECT_TESTS_DIR}/test_schurSolverPar.cpp
            #${PROJECT_TESTS_DIR}/test_schurNonOverlappingSolver.cpp
    )

    # Link test executable with the libraries to be tested and GoogleTest tests runner
    target_link_libraries(alltests PRIVATE
            NSBSolver
            gtest_main
    )

    # Include GoogleTest CMake module to enable test discovery
    include(GoogleTest)
    # Automatically discover and register all GoogleTest tests in the executable
    gtest_discover_tests(alltests)

    message(STATUS "Building tests - done")
endif ()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ========================================================
# ---- PROFILING -----------------------------------------
# ========================================================
# Define profiling-related variables

# Option to enable/disable profiling 
option(PROFILING "Build the tests" ON)
message(STATUS "PROFILING option: ${PROFILING}")
if (PROFILING)
    if (APPLE)
        # macOS: gprof non supportato. Abilita simboli per Instruments.
        message(STATUS "  -> macOS detected: Enabling Debug Symbols & Frame Pointer for Instruments.")
        set(PROFILING_FLAGS "-g" "-fno-omit-frame-pointer")
    else ()
        # Linux: Usa gprof standard
        message(STATUS "  -> Linux detected: Enabling '-pg' for gprof.")
        set(PROFILING_FLAGS "-pg")
    endif ()

    target_compile_options(NSBSolver PRIVATE ${PROFILING_FLAGS})
    target_link_libraries(NSBSolver PRIVATE ${PROFILING_FLAGS})

    target_compile_options(main PRIVATE ${PROFILING_FLAGS})
    target_link_libraries(main PRIVATE ${PROFILING_FLAGS})

    if (BUILD_TESTS)
        target_compile_options(alltests PRIVATE ${PROFILING_FLAGS})
        target_link_libraries(alltests PRIVATE ${PROFILING_FLAGS})
    endif ()


    message(STATUS "Building with profiling enabled - done")

    message(STATUS "Auto-profiling enabled: gmon.out -> profiling.txt")

    add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/.profile.stamp
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/main
            COMMAND gprof ${CMAKE_CURRENT_BINARY_DIR}/main gmon.out > ${CMAKE_CURRENT_BINARY_DIR}/profiling.txt
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            DEPENDS main
    )

    add_custom_target(profile
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/.profile.stamp
    )
endif ()