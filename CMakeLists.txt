cmake_minimum_required(VERSION 3.9)
project(Openbps C CXX)

set(CMAKE_CXX_STANDARD 17)

# Setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

#===============================================================================
# Command line options
#===============================================================================

option(optimize "Turn on O3 flag"        OFF)
option(openmp   "Enable shared-memory parallelism with OpenMP"   ON)
option(test     "Build test with google test package"   ON)

#===============================================================================
# Include source directories
#===============================================================================

include_directories(include/openbps)
include_directories(src)
include_directories(test)
include_directories(build)
include_directories(extern)

#===============================================================================
# Update git submodules as needed
#===============================================================================

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "git submodule update --init failed with \
        ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
    if(GIT_SUBMOD_RESULT EQUAL 0)
      message(STATUS "Eigen processing..")
      execute_process(COMMAND ${CMAKE_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")
    endif()
  endif()
endif()

add_definitions(-DHAVE_CBLAS=1)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)