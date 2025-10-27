cmake_minimum_required(VERSION 3.10)
project(NSBSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include headers
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
target_include_directories(InputReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/extern)

add_library(Initializer src/simulation/initializer.cpp)
target_include_directories(Initializer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Initializer PUBLIC Fields InputReader)

add_library(VTKWriter src/io/VTKWriter.cpp)
target_include_directories(VTKWriter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)


# Main executable
add_executable(test_initializer src/main.cpp)
target_link_libraries(test_initializer PRIVATE
    InputReader
    Initializer
    Fields
    VTKWriter
)