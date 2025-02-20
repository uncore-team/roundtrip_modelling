// filepath: /ransaconline2/ransaconline2/src/estimators/LN3Estimator.cpp
#include "LN3Estimator.h"
#include "Model.h"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cmath>

// LN3Estimator class implementation
// This class is responsible for fitting a log-normal model to the data
// and assessing its goodness of fit.

LN3Estimator::LN3Estimator() {
    // Constructor
}

LN3Estimator::~LN3Estimator() {
    // Destructor
}

Model LN3Estimator::fit(const std::vector<double>& data) {
    // Fit a log-normal model to the provided data.
    // DATA -> a vector of data to fit the model to.
    // Returns a Model object containing the fitted parameters.

    if (data.empty()) {
        throw std::invalid_argument("Data vector is empty.");
    }

    // Calculate the mean and standard deviation of the log-transformed data
    std::vector<double> logData(data.size());
    std::transform(data.begin(), data.end(), logData.begin(), [](double x) { return std::log(x); });

    double mean = std::accumulate(logData.begin(), logData.end(), 0.0) / logData.size();
    double variance = std::accumulate(logData.begin(), logData.end(), 0.0, [mean](double sum, double x) {
        return sum + (x - mean) * (x - mean);
    }) / logData.size();

    Model model;
    model.type = "LN3";
    model.defined = true;
    model.coeffs.mu = mean;
    model.coeffs.sigma = std::sqrt(variance);

    return model;
}

std::pair<bool, double> LN3Estimator::assess(const Model& model, const std::vector<double>& data) {
    // Assess the goodness of fit of the log-normal model to the data.
    // MODEL -> the fitted model to assess.
    // DATA -> the data to assess against the model.
    // Returns a pair containing a boolean indicating if the model fits well
    // and the statistic of the goodness of fit test.

    if (!model.defined) {
        throw std::invalid_argument("Model is not defined.");
    }

    // Implement the goodness of fit assessment logic here
    // For demonstration, we will return a dummy value
    double statistic = 0.0; // Replace with actual statistic calculation
    bool fitsWell = true; // Replace with actual fit assessment logic

    return {fitsWell, statistic};
}