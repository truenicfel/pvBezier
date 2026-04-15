# ------------------------------------
# --------- Project Setup ------------
# ------------------------------------

cmake_minimum_required(VERSION 3.18)

project(pvBezier
        LANGUAGES CXX
        DESCRIPTION "Reference Implementation for the EuroVis 2026 paper 'Parallel Vectors Extraction using Bézier Clipping'"
        HOMEPAGE_URL "https://vc.tf.fau.de/publications/Dassler26EuroVis/")

# ------------------------------------
# ------- Compiler Settings ----------
# ------------------------------------

set(CMAKE_CXX_STANDARD 23)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_BUILD_TYPE "Release")


# ------------------------------------
# ---------- Library Setup -----------
# ------------------------------------

include(FetchContent)

FetchContent_Declare(eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 9441d94dccccd5db8d64179516fdc5b53994a047
)
FetchContent_MakeAvailable(eigen3)

add_library(my_eigen INTERFACE)

target_link_libraries(my_eigen INTERFACE Eigen3::Eigen)

target_compile_definitions(my_eigen INTERFACE
        EIGEN_MATRIXBASE_PLUGIN="${CMAKE_CURRENT_SOURCE_DIR}/matrix_addons.hpp"
)

add_subdirectory(pvBezier)


# ------------------------------------
# ----------- Executable -------------
# ------------------------------------

add_executable(pvBezierMain main.cpp)
target_link_libraries(pvBezierMain pvBezierLib)
