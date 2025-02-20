cmake_minimum_required(VERSION 3.10)

# Project name and version
project(RansacOnline VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(ALG
    alglib/statistics.cpp
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
)

# Source files
set(SOURCES
    main.cpp
    Estimator.cpp
#    LL3Estimator.cpp
#   LN3Estimator.cpp
    ExponentialEstimator.cpp
)

# Include directories
include_directories(alglib)

## Add cmake target dependencies of the executable
## same as for the library above
## Specify libraries to link a library or executable target against

target_link_libraries(${PROJECT_NAME} ALG)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})