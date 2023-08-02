# CMake base setup
cmake_minimum_required(VERSION 3.19)
include(FetchContent)

# Define project info
project(
  NL-LatticeBoltzmann # change name here
  VERSION 0.0.1
  LANGUAGES CXX
  DESCRIPTION
  "lb - cdoe space"
  HOMEPAGE_URL "https://imtek-simulation.github.io/MolecularDynamics/"
)

# So that tests can work
enable_testing()

# We work with C++20 (-std=c++20)
set(CMAKE_CXX_STANDARD 20)
# trying out intrinic fkts
# march native enables everything already i guess lol
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -msse4.2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

# Try to find Eigen if it is installed
find_package(Eigen3 3.4.0 EXACT NO_MODULE)

# If not installed, fetch Eigen from the web
if(NOT Eigen3_FOUND)
  set(EIGEN_BUILD_TESTING OFF CACHE BOOL "Turn off tests" FORCE)

  FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
  )

  # We are not using FetchContent_MakeAvailable since that pulls in all
  # of Eigen's tests
  FetchContent_GetProperties(Eigen3)
  if (NOT Eigen3_POPULATED)
    FetchContent_Populate(Eigen3)  # Downloading Eigen
  endif()

  message(STATUS "Populated Eigen3 directory: ${eigen3_SOURCE_DIR}")

  # Setting up search of Eigen
  set(CMAKE_MODULE_PATH ${eigen3_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
  set(CMAKE_PREFIX_PATH ${eigen3_SOURCE_DIR} ${CMAKE_PREFIX_PATH})

  # Use FindEigen3.cmake from Eigen distribution
  find_package(Eigen3 3.4.0 REQUIRED MODULE EXACT)
endif()

# Option to turn ON/OFF MPI
option(USE_MPI "Activate MPI")

# Find MPI if requested
if(USE_MPI)
  find_package(MPI 3 REQUIRED)
endif()

# opencv
# Find OpenCV package
find_package(OpenCV REQUIRED)

# Include OpenCV headers
include_directories(${OpenCV_INCLUDE_DIRS})


add_subdirectory(simulation)  # contains the different iterations of the simulation code
add_subdirectory(src)         # Contains our common code pieces
add_subdirectory(tests)       # Tests for the library
add_subdirectory(milestones)  # Code for the different project milestones
