// filepath: /ransaconline2/ransaconline2/src/estimators/ExponentialEstimator.h
#ifndef EXPONENTIAL_ESTIMATOR_H
#define EXPONENTIAL_ESTIMATOR_H

#include <vector>
#include <optional>
#include <tuple>
#include "Estimator.h"

using namespace std;

// ExponentialEstimator class for fitting and assessing an exponential model to data.
class ExponentialEstimator : public Estimator {
public:
    // Constructor
    ExponentialEstimator();

    // Fit the exponential model to the provided data.
    // Parameters:
    // - data: A vector of double values representing the data to fit.
    // Returns: A boolean indicating success or failure of the fitting process.
    bool fit(const std::vector<double>& data, ModelParams&);

    // Assess the goodness of fit of the model to the data.
    // Parameters:
    // - data: A vector of double values representing the data to assess.
    // Returns: A tuple containing a boolean indicating if the model is rejected,
    //          the statistic of the goodness of fit test, and the threshold value.
    std::tuple<bool, double, double> assess(const std::vector<double>& data);

    // Get the fitted model.
    // Returns: An optional Model object containing the fitted model parameters.
    std::optional<Model> getModel() const;

private:
    int min_len;

};

#endif // EXPONENTIAL_ESTIMATOR_H