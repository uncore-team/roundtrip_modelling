#ifndef LOGNORMAL_ESTIMATOR_H
#define LOGNORMAL_ESTIMATOR_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Estimator.h"

using namespace std;

// LogNormalEstimator class for fitting and assessing an LogLogistic3 model to data.
class LogNormalEstimator : public Estimator {
public:
    // Constructor
    LogNormalEstimator();

    // Fit the LogLogistic model to the provided data.
    // Parameters:
    // - data: A vector of double values representing the data to fit.
    // Returns: A boolean indicating success or failure of the fitting process.
    Model fit(const vector<double>& samples) override;

    // Assess the goodness of fit of the model to the data.
    // Parameters:
    // - data: A vector of double values representing the data to assess.
    // Returns: A tuple containing a boolean indicating if the model is rejected,
    //          the statistic of the goodness of fit test, and the threshold value.
    tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) override;  // bool previous_model

    virtual double cdf(const ModelParams& params, const double& sample) override;
    virtual vector<double> cdf(const ModelParams& params, const vector<double>& samples) override;

    virtual double pdf(const ModelParams& params, const double& sample) override;
    virtual vector<double> pdf(const ModelParams& params, const vector<double>& samples) override;

    virtual double rnd(const ModelParams& params) override;
    virtual vector<double> rnd(const ModelParams& params, const unsigned& length) override;
};

#endif // LOGNORMAL_ESTIMATOR_H