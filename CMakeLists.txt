cmake_minimum_required(VERSION 3.10)

# Project name and version
project(OnlineRANSAC VERSION 1.0 LANGUAGES CXX)

# Enable multicore compilation
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N})
    message(STATUS "Using ${N} cores for parallel build")
endif()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force GCC compiler
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# Set default build type
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

####### Find required packages

# Find Python3
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
if(Python3_FOUND)
    # message(STATUS "Python3 found at ${Python3_EXECUTABLE}")
    include_directories(${Python3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Python3 not found. Please set CMAKE_PREFIX_PATH correctly.")
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    # message(STATUS "OpenMP found at ${OpenMP_gomp_LIBRARY}")
    include_directories(${OpenMP_INCLUDE_DIR})
else()
    message(FATAL_ERROR "OpenMP not found. Please set CMAKE_PREFIX_PATH correctly.")
endif()

# GTest library
file(GLOB GTest_SOURCES "${CMAKE_SOURCE_DIR}/gtest/src/*.cc")
add_library(GTest STATIC ${GTest_SOURCES})
target_include_directories(GTest PUBLIC "${CMAKE_SOURCE_DIR}/gtest" "${CMAKE_SOURCE_DIR}/gtest/include")

# ALGLIB library
file(GLOB ALGLIB_SOURCES "${CMAKE_SOURCE_DIR}/alglib/*.cpp")
add_library(ALGLIB STATIC ${ALGLIB_SOURCES})
target_include_directories(ALGLIB PUBLIC ${CMAKE_SOURCE_DIR}/alglib)

# CPU-specific optimization flags
if(MSVC OR MSYS OR MINGW)

    add_definitions(-DAE_OS=AE_WINDOWS)
    execute_process(COMMAND wmic CPU get Name OUTPUT_VARIABLE CPU_NAME)
    # message(STATUS "CPU detected: ${CPU_NAME}")

    if (${CPU_NAME} MATCHES "Intel")
        message(STATUS "Intel CPU detected, adding Intel-specific definitions")

		add_definitions(-DAE_CPU=AE_INTEL)
        if(MSVC) # Enable AVX2 for MSVC
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
        else() # Enable AVX2 for MinGW/MSYS
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
        endif()
    endif()

elseif(UNIX)

    add_definitions(-DAE_OS=AE_POSIX)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCHITECTURE)
    message(STATUS "Detected architecture: ${ARCHITECTURE}")

    if (${ARCHITECTURE} MATCHES "x86_64")
        message(STATUS "Intel CPU detected, adding Intel-specific definitions")
        add_definitions(-DAE_CPU=AE_INTEL)

        # Add AVX2 support flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
    endif()

endif()

# Main executable
set(SOURCES
    main.cpp
    Estimator.cpp
    OnlineRANSAC.cpp
    LogLogisticEstimator.cpp
    LogNormalEstimator.cpp
    ExponentialEstimator.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include
    ${CMAKE_CURRENT_SOURCE_DIR}/alglib
    ${Python3_INCLUDE_DIRS}
    ${OpenMP_INCLUDE_DIRS}
)

target_compile_options(${PROJECT_NAME} PRIVATE -libstdc++)
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Add compile definitions for proper I/O handling (needed for MinGW)
if(MINGW)
    add_definitions(
        -D_GLIBCXX_USE_CRT_STDIO_INPUT=1
        -D_GLIBCXX_USE_CRT_STDIO_OUTPUT=1
    )
endif()

target_link_options(${PROJECT_NAME} PRIVATE -static)
target_link_libraries(${PROJECT_NAME} PRIVATE 
    GTest
    ALGLIB
    gomp # OpenMP::OpenMP_CXX
    Python3::Python
    Python3::Module
    Python3::NumPy
    stdc++
)

# Testing configuration
enable_testing()

set(TEST_SOURCES
    tests/ExponentialTests.cpp
    # tests/LogLogisticTests.cpp
    tests/LogNormalTests.cpp
)

add_executable(${PROJECT_NAME}_tests 
    ${TEST_SOURCES}
    Estimator.cpp
    OnlineRANSAC.cpp
    LogLogisticEstimator.cpp
    LogNormalEstimator.cpp
    ExponentialEstimator.cpp
)

target_compile_options(${PROJECT_NAME}_tests PRIVATE -libstdc++)
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_link_options(${PROJECT_NAME}_tests PRIVATE -static)
target_link_libraries(${PROJECT_NAME}_tests PRIVATE 
    GTest
    ALGLIB
    gomp # OpenMP::OpenMP_CXX
    stdc++
)

target_include_directories(${PROJECT_NAME}_tests PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include
    ${CMAKE_CURRENT_SOURCE_DIR}/alglib
    ${OpenMP_INCLUDE_DIRS}
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}_tests)
