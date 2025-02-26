cmake_minimum_required(VERSION 3.10)

# Project name and version
project(OnlineRANSAC VERSION 1.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compile definitions for proper I/O handling (needed for MinGW)
add_definitions(
    -D_GLIBCXX_USE_CRT_STDIO_INPUT=1
    -D_GLIBCXX_USE_CRT_STDIO_OUTPUT=1
)

# Find OpenMP package
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    # Agregar flags de OpenMP a las opciones globales de compilación
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    
    # Para MinGW específicamente
    if(MINGW)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    endif()
endif()

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ALGLIB library sources
set(ALGLIB_SOURCES
    alglib/alglibinternal.cpp
    alglib/alglibmisc.cpp
    alglib/ap.cpp
    alglib/dataanalysis.cpp
    alglib/diffequations.cpp
    alglib/fasttransforms.cpp
    alglib/integration.cpp
    alglib/interpolation.cpp
    alglib/kernels_avx2.cpp
    alglib/kernels_fma.cpp
    alglib/kernels_sse2.cpp
    alglib/linalg.cpp
    alglib/optimization.cpp
    alglib/solvers.cpp
    alglib/specialfunctions.cpp
    alglib/statistics.cpp
)

# Create ALGLIB library
add_library(ALG STATIC ${ALGLIB_SOURCES})
target_include_directories(ALG PUBLIC ${CMAKE_SOURCE_DIR}/alglib)

# Main project sources
set(SOURCES
    main.cpp
    Estimator.cpp
    OnlineRANSAC.cpp
    LogLogisticEstimator.cpp
    LogNormalEstimator.cpp
    ExponentialEstimator.cpp
)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Configure include directories
target_include_directories(${PROJECT_NAME} 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/alglib
)

# Windows-specific configuration with MinGW
if(MINGW)
    target_link_options(${PROJECT_NAME} PRIVATE -static)
    target_link_libraries(${PROJECT_NAME} PRIVATE stdc++)
    target_compile_options(${PROJECT_NAME} PRIVATE -static-libstdc++ -static-libgcc)
endif()

# Configure warnings
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Link libraries
target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
        ALG
        OpenMP::OpenMP_CXX
)

# Configure OpenMP for ALGLIB
target_link_libraries(ALG PRIVATE OpenMP::OpenMP_CXX)
