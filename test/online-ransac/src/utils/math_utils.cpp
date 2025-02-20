#include "math_utils.h"
#include <cmath>

// Function to calculate the mean of a vector of doubles
double mean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}

// Function to calculate the variance of a vector of doubles
double variance(const std::vector<double>& data) {
    double data_mean = mean(data);
    double sum_squared_diff = 0.0;
    for (double value : data) {
        sum_squared_diff += (value - data_mean) * (value - data_mean);
    }
    return sum_squared_diff / data.size();
}

// Function to calculate the standard deviation of a vector of doubles
double standard_deviation(const std::vector<double>& data) {
    return std::sqrt(variance(data));
}

// Function to calculate the exponential probability density function
double exponential_pdf(double x, double lambda) {
    if (x < 0 || lambda <= 0) {
        return 0.0;
    }
    return lambda * std::exp(-lambda * x);
}