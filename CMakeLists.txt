cmake_minimum_required(VERSION 3.25)

project(project_1_fys4411)

# Set CXX standard to C++14
set(CMAKE_CXX_STANDARD 14)

# Check if the user has specified a build type via the command line, e.g.,
#
#   $ cmake .. -DCMAKE_BUILD_TYPE=Release
#
# if not, set Debug as the default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Add all compiler warnings
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# Add debugging flag
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add full compiler optimization when in the release build type
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# ==== CHANGE 1: Find OpenMP package ====
find_package(OpenMP REQUIRED)

# ==== CHANGE 2: Find MPI package ====
find_package(MPI REQUIRED)

# ==== CHANGE 3: Find Eigen package ====
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Add directories containing header-files
include_directories(
    .
    Hamiltonians
    InitialStates
    Math
    WaveFunctions
    Solvers
)

# Add sources from all directories
file(
    GLOB SOURCES "*.cpp"
    GLOB SOURCES "Hamiltonians/*.cpp"
    GLOB SOURCES "InitialStates/*.cpp"
    GLOB SOURCES "Math/*.cpp"
    GLOB SOURCES "Solvers/*.cpp"
    GLOB SOURCES "WaveFunctions/*.cpp"
)

# Create executable called "vmc"
add_executable(vmc ${SOURCES})

# ==== CHANGE 4: Link OpenMP to the executable ====
if(OpenMP_CXX_FOUND)
    target_link_libraries(vmc PUBLIC OpenMP::OpenMP_CXX)
endif()

# ==== CHANGE 5: Link MPI to the executable ====
if(MPI_CXX_FOUND)
    target_link_libraries(vmc PUBLIC MPI::MPI_CXX)
endif()
