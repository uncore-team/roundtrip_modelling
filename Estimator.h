#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <memory>
#include <vector>
#include <tuple>

#include "Model.h"

using namespace std;

// Estimator abstract class for fitting some model to data
class Estimator {
public:
    using Ptr = shared_ptr<Estimator>;

    Estimator(); // Constructor
    virtual ~Estimator();  // Destructor virtual

    // Fit the log-normal model to the provided data
    // Parameters:
    // - data: A vector of double values representing the data to fit
    // Returns: A Model object containing the fitted model parameters
    virtual Model fit(const vector<double>& samples) = 0;

    // Assess the goodness of fit for the fitted model
    // Parameters:
    // - data: A vector of double values representing the data to assess
    // Returns: A tuple containing a boolean indicating rejection of the null hypothesis,
    //          the statistic value, and the threshold for the goodness of fit test
    virtual std::tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) = 0;
};

#endif // ESTIMATOR_H