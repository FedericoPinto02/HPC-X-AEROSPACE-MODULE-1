cmake_minimum_required(VERSION 3.15)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ========================================================
# ---- SOLVER LIBRARY ------------------------------------
# ========================================================
# Define project-related variables
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(PROJECT_EXTERN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern) # for nlohmann/json

# Fetch MuParser source code automatically
FetchContent_Declare(
        muparserx
        GIT_REPOSITORY https://github.com/beltoforion/muparserx.git
        GIT_TAG v4.0.12
)
FetchContent_MakeAvailable(muparserx)

file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/*.cpp)
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*main.cpp$") # exclude the main.cpp
file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS ${PROJECT_INCLUDE_DIR}/*.hpp)

# Add owned src and headers in the library
add_library(NSBSolver ${SOURCE_FILES} ${HEADER_FILES})
# Include owned headers (exposed publicly), external headers (internal use only)
# and muparserx headers (internal use only).
target_include_directories(NSBSolver
        PUBLIC ${PROJECT_INCLUDE_DIR}
        PRIVATE ${PROJECT_EXTERN_DIR} ${muparserx_SOURCE_DIR}/parser
)
# Link muparserx library for internal use only
target_link_libraries(NSBSolver PRIVATE muparserx)


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
            ${PROJECT_TESTS_DIR}/test_tridiag.cpp
            ${PROJECT_TESTS_DIR}/test_fields.cpp
            ${PROJECT_TESTS_DIR}/test_derivatives.cpp
            ${PROJECT_TESTS_DIR}/test_linearSys.cpp
            ${PROJECT_TESTS_DIR}/test_viscousStep.cpp
            ${PROJECT_TESTS_DIR}/test_pressureStep.cpp
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