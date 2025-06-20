cmake_minimum_required(VERSION 3.15)
project(meshgen LANGUAGES CXX)

# -----------------------------
# Compiler and language settings
# -----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable OpenMP and parallel STL if using GCC
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Using GCC: enabling parallel STL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -D_GLIBCXX_PARALLEL")
endif()

# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories("/home/debasish/Desktop/softwares/eigen-3.4.0")

# ------------------------
# OpenMP support
# ------------------------
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP")
endif()

# ------------------------
# Main executable: meshgen
# ------------------------
set(SRC_FILES
    main/main.cpp
    src/write_txt.cpp
    src/write_vtk.cpp
    src/mesh_utils.cpp
    src/quadrature.cpp
    src/deduplicate.cpp
    src/connectivity.cpp
    src/geometry_analyzer.cpp
    src/geometry_utility.cpp
    src/trgl6_icos.cpp
)

add_executable(meshgen ${SRC_FILES})
target_link_libraries(meshgen PRIVATE OpenMP::OpenMP_CXX)

# ------------------------
# Final summary
# ------------------------
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Output binary directory: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
