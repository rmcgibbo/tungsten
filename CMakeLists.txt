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
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
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
add_executable(tungsten ${TUNGSTEN_SOURCES})

# Testsing
enable_testing()
add_library(libtungsten ${TUNGSTEN_SOURCES})
add_executable(TestReductions tests/TestReductions.cpp)
target_link_libraries(TestReductions libtungsten)
set(TARGETS tungsten TestReductions)


# Our includes
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/rmsd/")
# MPI includes
include_directories(${MPI_CXX_INCLUDE_PATH})
# OpenMM includes
include_directories(${OpenMM_INCLUDE_PATH})
# NetCDF includes
include_directories(${NetCDF_INCLUDE_DIR})
include_directories(${NetCDF_INCLUDE_DIRS})

foreach(target ${TARGETS})
    # MPI link
    target_link_libraries(${target} ${MPI_CXX_LIBRARIES})
    target_link_libraries(${target} ${MPI_CXX_LIBRARIES})
    # OpenMM /link
    target_link_libraries(${target} ${OpenMM_LIBRARY})
    target_link_libraries(${target} ${OpenMM_SERIALIZATION_LIBRARY})
    # NetCDF link
    target_link_libraries(${target} ${NetCDF_C_LIBRARY})
endforeach(target ${TARGETS})

add_test(TestReductions2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 TestReductions)
add_test(TestReductions5 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 TestReductions)

