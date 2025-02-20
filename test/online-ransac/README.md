# online-ransac/README.md

# Online RANSAC

## Overview

The Online RANSAC project implements an online version of the RANSAC (Random Sample Consensus) algorithm specifically designed for estimating exponential distributions. This project provides a robust framework for modeling and updating parameters in real-time as new data becomes available.

## Features

- **Online RANSAC Algorithm**: Efficiently estimates parameters of exponential distributions using the RANSAC approach.
- **Model Fitting**: Automatically fits the exponential model to incoming data.
- **Parameter Updating**: Dynamically updates model parameters as new data points are processed.

## Installation

To build the project, ensure you have CMake installed. Clone the repository and navigate to the project directory:

```bash
git clone <repository-url>
cd online-ransac
```

Then, create a build directory and compile the project:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

After building the project, you can run the application using the following command:

```bash
./online-ransac
```

You can modify the `src/main.cpp` file to customize the behavior of the Online RANSAC algorithm according to your needs.

## Testing

Unit tests are provided for both the `OnlineRANSAC` and `ExponentialModel` classes. To run the tests, navigate to the build directory and execute:

```bash
make test
```

## Contributing

Contributions are welcome! Please submit a pull request or open an issue for any enhancements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.