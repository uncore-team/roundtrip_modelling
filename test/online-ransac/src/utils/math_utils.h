#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>

// Function to calculate the mean of a vector of doubles
double calculateMean(const std::vector<double>& data);

// Function to calculate the variance of a vector of doubles
double calculateVariance(const std::vector<double>& data);

// Function to calculate the standard deviation of a vector of doubles
double calculateStandardDeviation(const std::vector<double>& data);

#endif // MATH_UTILS_H