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
      execute_process(COMMAND cmake ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen)
    endif()
  endif()
endif()

#===============================================================================
# Include extern packages
#===============================================================================

add_subdirectory(extern/xtl)
add_subdirectory(extern/xtensor)
target_link_libraries(xtensor INTERFACE xtl)
add_subdirectory(lib/googletest)
add_definitions(-DHAVE_CBLAS=1)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#Binary build
#===============================================================================
# Compiler options
#===============================================================================
#Check compiler options
if(optimize)
  list(APPEND cxxflags -O3)
endif()
if(openmp)
  # Requires CMake 3.1+
  find_package(OpenMP)
  if(OPENMP_FOUND)
    list(APPEND cxxflags ${OpenMP_CXX_FLAGS})
    list(APPEND ldflags ${OpenMP_CXX_FLAGS})
  endif()
endif()
# Show flags being used
message(STATUS "C++ flags: ${cxxflags}")

#===============================================================================
# Setup process
#===============================================================================
#Linking with thread library
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -lrt")

list(APPEND lib_SOURCES
      "src/chain.cpp"
      "src/configure.cpp"
      "src/executer.cpp"
      "src/filter.cpp"
      "src/functionals.cpp"
      "src/materials.cpp"
      "src/matrix.cpp"
      "src/nuclide.cpp"
      "src/parse.cpp"
      "src/reactions.cpp"
      "extern/pugiData/pugixml.cpp")

add_library(openbps_lib STATIC ${lib_SOURCES})
# Includes
target_include_directories(openbps_lib
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
# Set compile flags
target_compile_options(openbps_lib PRIVATE ${cxxflags})
# Link library
target_link_libraries(openbps_lib Eigen3::Eigen xtensor)

add_executable(openbps "../src/main.cpp" )
# Set compile flags
target_compile_options(openbps PRIVATE ${cxxflags})
# Link executable
target_link_libraries(openbps PUBLIC openbps_lib)

#Test build
if(test)
add_executable(o_tst  "../test/test.cpp"
                      "../test/main_test.cpp"
                      "../lib/googletest/googletest/include/gtest/gtest.h"
                      ${lib_SOURCES})

add_test(NAME o_tst COMMAND o_tst)
target_link_libraries(o_tst PUBLIC openbps_lib gtest)
endif()

