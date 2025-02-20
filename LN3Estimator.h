// filepath: /ransaconline2/ransaconline2/src/estimators/LN3Estimator.h
#ifndef LN3ESTIMATOR_H
#define LN3ESTIMATOR_H

#include <vector>
#include <optional>
#include "Model.h"

// LN3Estimator class for fitting a log-normal model to data
class LN3Estimator {
public:
    // Constructor
    LN3Estimator();

    // Fit the log-normal model to the provided data
    // Parameters:
    // - data: A vector of double values representing the data to fit
    // Returns: A Model object containing the fitted model parameters
    Model fit(const std::vector<double>& data);

    // Assess the goodness of fit for the fitted model
    // Parameters:
    // - data: A vector of double values representing the data to assess
    // Returns: A tuple containing a boolean indicating rejection of the null hypothesis,
    //          the statistic value, and the threshold for the goodness of fit test
    std::tuple<bool, double, double> assess(const std::vector<double>& data);

private:
    Model model; // The fitted log-normal model
};

#endif // LN3ESTIMATOR_H