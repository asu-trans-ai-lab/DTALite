cmake_minimum_required (VERSION 3.8)

project ("DTALite")

option(BUILD_EXE "BUILD Executable or not" ON)

configure_file(build_config.h.in build_config.h)

if(BUILD_EXE)
    message("Building executable on ${CMAKE_SYSTEM_NAME}")
else()
    message("Building shared library on ${CMAKE_SYSTEM_NAME}")
endif()

# specify the C++ standard
# C++11 is necessary as teestream.h requires C++11 to be complied
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# enable OpenMP
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp -fpic")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-fopenmp -fpic")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # do nothing as MSVC does not support -fopenmp -fpic
else()
    message("Please set OpenMP flag compatible with your compiler")
endif()

# note CMAKE_BUILD_TYPE WILL BE INGNORED for IDE generators, like Visual Studio and Xcode
# it only works for single-configuration generators, e.g., make and Ninja
# see https://cmake.org/cmake/help/v3.0/variable/CMAKE_BUILD_TYPE.html?highlight=cmake_build_types for details
set(CMAKE_BUILD_TYPE  "Release")
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/bin)

if(BUILD_EXE)
    add_executable(DTALite utils.cpp flash_dta.cpp main_api.cpp simulation.cpp)
else()
    add_library(DTALite SHARED utils.cpp main_api.cpp simulation.cpp)
endif()

# using OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(DTALite PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(DTALite PUBLIC "${PROJECT_BINARY_DIR}")