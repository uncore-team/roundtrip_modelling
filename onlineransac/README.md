# ransaconline# OnlineRANSAC

Implementation of an Online RANSAC algorithm with different estimators.

## Project Structure

```
.
├── CMakeLists.txt          # Main build configuration file
├── CMakePresets.json       # CMake presets configuration
├── src/                    # Source files
│   ├── Estimator.*        # Base estimator class
│   ├── ExponentialEstimator.*   # Exponential distribution estimator
│   ├── LogLogisticEstimator.*   # Log-logistic distribution estimator
│   ├── LogNormalEstimator.*     # Log-normal distribution estimator
│   ├── OnlineRANSAC.*     # Main RANSAC implementation
│   ├── Model.h            # Model interface
│   ├── misc.cpp           # Utility functions
│   └── matplotlibcpp.h    # Plotting library header
├── tests/                  # Test directory
├── alglib/                 # Algebraic algorithms library
├── gtest/                  # Google Test framework
├── matlab/                 # MATLAB scripts and implementations
└── .vscode/               # VS Code configuration files
```

## Dependencies

- CMake (build system)
- Matplotlib for C++ (plotting)
- Google Test (testing framework)
- ALGLIB (numerical analysis and data processing)

## Building the Project

The project uses CMake as its build system. See CMakeLists.txt for build configuration details.
