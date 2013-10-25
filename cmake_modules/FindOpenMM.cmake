# Find OpenMM library.
#
# Looks for the OpenMM libraries at the default (/usr/local) location 
# or custom location found in the OPENMM_DIR, OPENMM_INCLUDE_PATH
# or OPENMM_LIB_PATH environment variables. 
#
# The script defines defines: 
#  OpenMM_FOUND     
#  OpenMM_DIR
#  OpenMM_INCLUDE_PATH
#  OpenMM_LIB_PATH
#  OpenMM_LIBRARY
#  OpenMM_SERIALIZATION_LIBRARY      
#  OpenMM_PLUGIN_DIR
#  OpenMM_HAVE_VERSION_52_PLUS
#
# Author: Szilard Pall (pszilard@cbr.su.se)
# Modified by Robert McGibbon (rmcgibbo@gmail.com)

include(PrintVariable)
include(CheckIncludeFileCXX)
include(FindPackageHandleStandardArgs)

if(OpenMM_INCLUDE_PATH AND OpenMM_LIB_PATH AND OpenMM_PLUGIN_DIR)
    set(OpenMM_FIND_QUIETLY)
endif()

# try getting OpenMM root directory from environment variable or default
# location 
file(TO_CMAKE_PATH "$ENV{OPENMM_DIR}" _env_OPENMM_DIR)
if(IS_DIRECTORY ${_env_OPENMM_DIR})
    set(OpenMM_DIR "${_env_OPENMM_DIR}")
else()
    if(IS_DIRECTORY "/usr/local/openmm")
        set(OpenMM_DIR "/usr/local/openmm")
    endif()
endif()

# try looking for the environment variable OPENMM_INCLUDE_PATH, which the
# python installer for openmm requires, so it might be around
if (NOT DEFINED OpenMM_DIR)
  file(TO_CMAKE_PATH "$ENV{OPENMM_INCLUDE_PATH}" _env_OPENMM_INCLUDE_PATH)
	if(IS_DIRECTORY ${_env_OPENMM_INCLUDE_PATH})
	    #set(OpenMM_INCLUDE_PATH  "${_env_OPENMM_INCLUDE_PATH}")
	    set(OpenMM_DIR "${_env_OPENMM_INCLUDE_PATH}/..")
	endif()
endif()

# try looking for the OPENMM_LIB_PATH environment variable, which the
# python installer also might use
if (NOT DEFINED OpenMM_DIR)
  file(TO_CMAKE_PATH "$ENV{OPENMM_LIB_PATH}" _env_OPENMM_LIB_PATH)
	if(IS_DIRECTORY ${_env_OPENMM_LIB_PATH})
	    #set(OpenMM_LIB_PATH  "${_env_OPENMM_LIB_PATH}")
	    set(OpenMM_DIR "${_env_OPENMM_LIB_PATH}/..")
	endif()
endif()

# set all of the paths, based on whichever way we found them
if (DEFINED OpenMM_DIR)
	get_filename_component(OpenMM_DIR  "${OpenMM_DIR}" REALPATH)
	get_filename_component(OpenMM_LIB_PATH  "${OpenMM_DIR}/lib" REALPATH)
	get_filename_component(OpenMM_INCLUDE_PATH  "${OpenMM_DIR}/include" REALPATH)
endif()

find_library(OpenMM_LIBRARY
    NAMES OpenMM
    PATHS "${OpenMM_LIB_PATH}"
    CACHE STRING "OpenMM library")

find_library(OpenMM_SERIALIZATION_LIBRARY
    NAMES OpenMMSerialization
    PATHS "${OpenMM_LIB_PATH}"
    CACHE STRING "OpenMM serialization library")


find_path(OpenMM_INCLUDE_PATH 
    NAMES OpenMM.h 
    PATHS "${OpenMM_DIR}/include" "${OpenMM_LIB_PATH}/../include" "${OpenMM_INCLUDE_PATH}")


if(NOT IS_DIRECTORY ${OpenMM_DIR})
    message(FATAL_ERROR "Could not find OpenMM! Set the OPENMM_DIR environment "
    "variable to contain the path of the OpenMM installation.")
endif()

if(NOT IS_DIRECTORY ${OpenMM_LIB_PATH})
    message(FATAL_ERROR "Can't find OpenMM libraries. Check your OpenMM installation!")
endif()

if(NOT OpenMM_SERIALIZATION_LIBRARY)
   message(ERROR "Could not find OpenMM serialization library. Tungsten requires that OpenMM be built with serialization support.")
endif()

if(NOT OpenMM_INCLUDE_PATH)
    message(FATAL_ERROR "Can't find OpenMM includes. Check your OpenMM installation!")
endif()

if(IS_DIRECTORY "${OpenMM_LIB_PATH}/plugins")
    get_filename_component(OpenMM_PLUGIN_DIR
        "${OpenMM_LIB_PATH}/plugins"
        ABSOLUTE)
    set(OpenMM_PLUGIN_DIR ${OpenMM_PLUGIN_DIR} CACHE PATH "OpenMM plugins directory")
else()
    message(WARNING "Could not detect the OpenMM plugin directory at the default location (${OpenMM_LIB_PATH}/plugins). Check your OpenMM installation.")
endif()



set(CMAKE_REQUIRED_INCLUDES ${OpenMM_INCLUDE_PATH} )
check_include_file_cxx(openmm/MonteCarloAnisotropicBarostat.h
                       OpenMM_HAVE_VERSION_52_PLUS)
set(CMAKE_REQUIRED_INCLUDES)

set(OpenMM_DIR ${OpenMM_DIR} CACHE PATH "OpenMM installation directory")
set(OpenMM_PLUGIN_DIR ${OpenMM_PLUGIN_DIR} CACHE PATH "OpenMM plugin directory")
set(OpenMM_INCLUDE_PATH ${OpenMM_INCLUDE_PATH} CACHE PATH "OpenMM include directory")



find_package_handle_standard_args(OpenMM DEFAULT_MSG 
                                    OpenMM_DIR
                                    OpenMM_LIBRARY
                                    OpenMM_SERIALIZATION_LIBRARY
                                    OpenMM_LIB_PATH 
                                    OpenMM_INCLUDE_PATH)


mark_as_advanced(
  OpenMM_PLUGIN_DIR
)