cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 OLD)  # the netcdf find package seems to require this
project (tungsten)

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
find_package(MPI)
find_package(NetCDF)
find_package(OpenMM)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Platform Cheks
include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckTypeSize)
check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
check_function_exists(gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists(_aligned_malloc HAVE_ALIGNED_MALLOC)
check_include_file_cxx(sys/time.h HAVE_SYS_TIME_H)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/config.h)

# set the source directory
file(GLOB_RECURSE TUNGSTEN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*")
add_executable (tungsten ${TUNGSTEN_SOURCES})

# Our includes
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/rmsd/")

# MPI include/link
target_link_libraries(tungsten ${MPI_CXX_LIBRARIES})
target_link_libraries(tungsten ${MPI_CXX_LIBRARIES})
include_directories(${MPI_CXX_INCLUDE_PATH})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# OpenMM include/link
include_directories(${OpenMM_INCLUDE_PATH})
target_link_libraries(tungsten ${OpenMM_LIBRARY})
target_link_libraries(tungsten ${OpenMM_SERIALIZATION_LIBRARY})

# NetCDF include/link
include_directories(${NetCDF_INCLUDE_DIR})
include_directories(${NetCDF_INCLUDE_DIRS})
target_link_libraries(tungsten ${NetCDF_C_LIBRARY})

