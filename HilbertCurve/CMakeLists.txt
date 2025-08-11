cmake_minimum_required(VERSION 3.19.2)

project(TorusToHilbertCurve LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin) # For potential future static/shared libs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin) # For potential future static libs

# Include directories (Global - consider using target_include_directories more)
# include_directories(${CMAKE_SOURCE_DIR}/include) # Can be removed if all libs use target_include_directories

# Define libraries
# Torus library (header-only)
add_library(torus INTERFACE)
target_include_directories(torus INTERFACE ${CMAKE_SOURCE_DIR}/include)

# Torus Donut library (header-only)
add_library(torus_donut INTERFACE)
target_include_directories(torus_donut INTERFACE ${CMAKE_SOURCE_DIR}/include)

# Hilbert Curve library (now header-only)
add_library(hilbert_curve INTERFACE)
target_include_directories(hilbert_curve INTERFACE ${CMAKE_SOURCE_DIR}/include)

# Mesh library (now header-only)
add_library(mesh INTERFACE)
target_include_directories(mesh INTERFACE ${CMAKE_SOURCE_DIR}/include)

# Dist library (header-only)
add_library(dist INTERFACE)
target_include_directories(dist INTERFACE ${CMAKE_SOURCE_DIR}/include)

# Define the main executable
add_executable(Torus_HilbertCurve ${CMAKE_SOURCE_DIR}/src/main.cpp)

# Link libraries to the executable
# Linking INTERFACE libraries correctly propagates include directories
target_link_libraries(Torus_HilbertCurve PRIVATE torus torus_donut hilbert_curve mesh dist)

# Define a cache variable for runtime arguments
set(RUN_ARGS "" CACHE STRING "Arguments to pass to Torus_HilbertCurve when using run_all")

# Custom target to build and run the executable with args
add_custom_target(run_all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target Torus_HilbertCurve
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Torus_HilbertCurve ${RUN_ARGS}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Building and running Torus_HilbertCurve with arguments: ${RUN_ARGS}"
    DEPENDS Torus_HilbertCurve # Ensure executable is built first
)

# Optional: Add installation rules
install(TARGETS Torus_HilbertCurve
    RUNTIME DESTINATION bin
)
# Install all headers from the include directory
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)