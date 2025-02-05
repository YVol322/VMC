cmake_minimum_required(VERSION 3.25)

project(project_1_fys4411)

# Set CXX standard to C++14
set(CMAKE_CXX_STANDARD 14)

# Check if the user has specified a build type via the command line, e.g., (in the
# build-directory)
#
#   $cmake .. -DCMAKE_BUILD_TYPE=Release
#
# if not, set Debug as the default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Add all compiler warnings. Remember, a warning in C++ can very often lead to a bug
# later on so _fix your warnings!_
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# Add debugging flag for tools such as gdb and valgrind.
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add full compiler optimization when in the release build type. Other options include
# -O0, -O1, -O2.
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# ==== CHANGE 1: Find OpenMP package ====
find_package(OpenMP REQUIRED)

# ==== CHANGE 2: Find MPI package ====
find_package(MPI REQUIRED)

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

# ==== CHANGE 3: Link OpenMP to the executable ====
if(OpenMP_CXX_FOUND)
    target_link_libraries(vmc PUBLIC OpenMP::OpenMP_CXX)
endif()

# ==== CHANGE 4: Link MPI to the executable ====
if(MPI_CXX_FOUND)
    target_link_libraries(vmc PUBLIC MPI::MPI_CXX)
endif()