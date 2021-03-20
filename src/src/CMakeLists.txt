cmake_minimum_required (VERSION 3.8)

project ("DTALite")

option(BUILD_EXE "BUILD Executable or not" ON)

configure_file(config.h.in config.h)

if(BUILD_EXE)
    message("Building executable on ${CMAKE_SYSTEM_NAME}")
elseif()
    message("Building shared library on ${CMAKE_SYSTEM_NAME}")
endif()

# specify the C++ standard
# C++11 is necessary as teestream.h requires C++11 to be complied
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# enable OpenMP
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp -lomp fpic")
elseif()
    set(CMAKE_CXX_FLAGS "-fopenmp -fpic")
endif()

# note CMAKE_BUILD_TYPE WILL BE INGNORED for IDE generators, like Visual Studio and Xcode
# it only works for single-configuration generators, e.g., make and Ninja
# see https://cmake.org/cmake/help/v3.0/variable/CMAKE_BUILD_TYPE.html?highlight=cmake_build_types for details
set(CMAKE_BUILD_TYPE  "Release")
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/bin)

if(BUILD_EXE)
    add_executable(DTALite flash_dta.cpp main_api.cpp)
elseif()
    add_library(DTALite SHARED main_api.cpp)
endif()

# using OpenMP, the second command is only needed for macOS
find_package(OpenMP REQUIRED)
target_link_libraries(DTALite PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(DTALite PUBLIC "${PROJECT_BINARY_DIR}")