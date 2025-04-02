cmake_minimum_required(VERSION 3.10)
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

# # Find Python3 (dynamic linking)
# find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
# if(Python3_FOUND)
#     message(STATUS "Python3 found.")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_PYTHON")
# else()
#     message(WARNING "Python3 not found. Please set CMAKE_PREFIX_PATH correctly.")
# endif()

# Find OpenMP (dynamic linking)
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "OpenMP found.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP not found. Please set CMAKE_PREFIX_PATH correctly.")
endif()

# GTest (static linking)
file(GLOB GTest_SOURCES "${CMAKE_SOURCE_DIR}/gtest/src/*.cc")
add_library(GTest STATIC ${GTest_SOURCES})
target_include_directories(GTest PUBLIC "${CMAKE_SOURCE_DIR}/gtest" "${CMAKE_SOURCE_DIR}/gtest/include")

# ALGLIB (static linking)
file(GLOB ALGLIB_SOURCES "${CMAKE_SOURCE_DIR}/alglib/*.cpp")
add_library(ALGLIB STATIC ${ALGLIB_SOURCES})
target_include_directories(ALGLIB PUBLIC ${CMAKE_SOURCE_DIR}/alglib)

# CPU-specific optimization flags
if(MSVC OR MSYS OR MINGW)
    add_definitions(-DAE_OS=AE_WINDOWS)
    execute_process(COMMAND wmic CPU get Name OUTPUT_VARIABLE CPU_NAME)
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
    if (${ARCHITECTURE} MATCHES "x86_64")
        message(STATUS "Intel CPU detected, adding Intel-specific definitions")
        add_definitions(-DAE_CPU=AE_INTEL)
        # Add AVX2 support flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
    endif()
endif()

# Main executable
set(SOURCES
    misc.cpp
    Estimator.cpp
    OnlineRANSAC.cpp
    LogLogisticEstimator.cpp
    LogNormalEstimator.cpp
    ExponentialEstimator.cpp
)

add_executable(${PROJECT_NAME}
    main.cpp
    ${SOURCES}
)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include
    ${CMAKE_CURRENT_SOURCE_DIR}/alglib
    # ${Python3_INCLUDE_DIRS}
    ${OpenMP_INCLUDE_DIRS}
)

# Compiler options
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic -Wunused)
endif()

# Add compile definitions for proper I/O handling (needed for MinGW)
if(MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

# Link libraries (static and dynamic)
target_link_libraries(${PROJECT_NAME} PRIVATE 
    # Static
    GTest
    ALGLIB
    # Dynamic
    # Python3::Python
    # Python3::Module
    # Python3::NumPy
    OpenMP::OpenMP_CXX
)

# Testing configuration
enable_testing()

set(TEST_SOURCES
    tests/ExponentialTests.cpp
    tests/LogLogisticTests.cpp
    tests/LogNormalTests.cpp
)

# Test executable
add_executable(${PROJECT_NAME}_tests 
    ${TEST_SOURCES}
    ${SOURCES}
)

# Include directories for tests
target_include_directories(${PROJECT_NAME}_tests PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/gtest/include
    ${CMAKE_CURRENT_SOURCE_DIR}/alglib
    ${OpenMP_INCLUDE_DIRS}
)

# Compiler options for tests
if(MSVC)
    target_compile_options(${PROJECT_NAME}_tests PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME}_tests PRIVATE -Wall -Wextra -Wpedantic -Wunused)
endif()

# Add compile definitions for proper I/O handling (needed for MinGW)
if(MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

# Link libraries for tests (static and dynamic)
target_link_libraries(${PROJECT_NAME}_tests PRIVATE 
    # Static
    GTest
    ALGLIB
    # Dynamic
    OpenMP::OpenMP_CXX
)

# GoogleTest integration
include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}_tests)
