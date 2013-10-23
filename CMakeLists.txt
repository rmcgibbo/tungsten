cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 OLD)  # the netcdf find package seems to require this

project (tungsten)

# Set the compilers to mpicc
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Custom cmake modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

# set the default build type as release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Find all of our dependencies
find_package(NetCDF)
find_package(OpenMM)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


file(GLOB_RECURSE TUNGSTEN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*")
add_executable (tungsten ${TUNGSTEN_SOURCES})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/rmsd/")
include_directories(${OpenMM_INCLUDE_PATH})
include_directories(${NetCDF_INCLUDE_DIR})
include_directories(${NetCDF_INCLUDE_DIRS})

target_link_libraries(tungsten ${NetCDF_CXX_LIBRARY})
target_link_libraries(tungsten ${OpenMM_LIBRARY})
target_link_libraries(tungsten ${OpenMM_SERIALIZATION_LIBRARY})

