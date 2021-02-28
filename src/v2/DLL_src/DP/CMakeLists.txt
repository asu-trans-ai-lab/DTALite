cmake_minimum_required (VERSION 3.8)

project ("DTALite_DLL")

if(WIN32)
    message("Building dll file on Windows")
elseif(UNIX)
    message("Building so or dylib file on ${CMAKE_SYSTEM_NAME}")
endif()

# enable OpenMP
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp -lomp")
elseif()
    set(CMAKE_CXX_FLAGS "-fopenmp -fpic")
endif()

# note CMAKE_BUILD_TYPE WILL BE INGNORED for IDE generators, like Visual Studio and Xcode
# it only works for single-configuration generators, e.g., make and Ninja
# see https://cmake.org/cmake/help/v3.0/variable/CMAKE_BUILD_TYPE.html?highlight=cmake_build_types for details
set(CMAKE_BUILD_TYPE  "Release")
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/bin)

add_library(DTALite_DLL SHARED MathLibrary.cpp main_api.cpp)

# using OpenMP, the second command is only needed for macOS
find_package(OpenMP REQUIRED)
target_link_libraries(DTALite_DLL PRIVATE OpenMP::OpenMP_CXX)