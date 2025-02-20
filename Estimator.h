#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <vector>
#include "Model.h"

using namespace std;

// LN3Estimator class for fitting a log-normal model to data
class Estimator {
public:
    // Constructor
    Estimator() {};

    // Fit the log-normal model to the provided data
    // Parameters:
    // - data: A vector of double values representing the data to fit
    // Returns: A Model object containing the fitted model parameters
    virtual Model fit(const vector<double>& data) = 0;

    // Assess the goodness of fit for the fitted model
    // Parameters:
    // - data: A vector of double values representing the data to assess
    // Returns: A tuple containing a boolean indicating rejection of the null hypothesis,
    //          the statistic value, and the threshold for the goodness of fit test
    virtual std::tuple<bool, double, double> assess(const vector<double>& data) = 0;
};

#endif // ESTIMATOR_H