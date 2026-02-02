cmake_minimum_required(VERSION 3.16)
project(magcal_cpp_full LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(MAGCAL_ENABLE_MATLAB_ENGINE "Enable MATLAB Engine wrapper for SOS (requires MATLAB)" OFF)

find_package(Eigen3 3.3 QUIET)

if(NOT Eigen3_FOUND)
    # Fallback to a simple header search. Many systems install Eigen headers at /usr/include/eigen3.
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core PATH_SUFFIXES eigen3)
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR "Eigen3 was not found. Install Eigen (e.g., libeigen3-dev on Ubuntu) or set Eigen3_DIR / EIGEN3_INCLUDE_DIR.")
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

add_library(magcal
    src/magcal.cpp
)

target_include_directories(magcal PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(magcal PUBLIC Eigen3::Eigen)

target_compile_definitions(magcal PUBLIC EIGEN_MPL2_ONLY)

target_compile_options(magcal PRIVATE -Wall -Wextra -Wpedantic)

add_executable(magcal_run src/main.cpp)

target_link_libraries(magcal_run PRIVATE magcal)

if(MAGCAL_ENABLE_MATLAB_ENGINE)
    message(STATUS "MATLAB Engine wrapper requested. You must configure MATLAB include/lib paths manually.")
    target_compile_definitions(magcal PUBLIC MAGCAL_ENABLE_MATLAB_ENGINE)
endif()
