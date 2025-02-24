#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <memory>
#include <vector>
#include <tuple>

#include "Model.h"

using namespace std;

// LN3Estimator class for fitting a log-normal model to data
class Estimator {
public:
    using Ptr = std::shared_ptr<Estimator>;

    // Constructor
    Estimator();
    virtual ~Estimator() = default;  // Destructor virtual

    // Fit the log-normal model to the provided data
    // Parameters:
    // - data: A vector of double values representing the data to fit
    // Returns: A Model object containing the fitted model parameters
    virtual Model fit(const std::vector<double>& samples) = 0;  // Método virtual

    // Assess the goodness of fit for the fitted model
    // Parameters:
    // - data: A vector of double values representing the data to assess
    // Returns: A tuple containing a boolean indicating rejection of the null hypothesis,
    //          the statistic value, and the threshold for the goodness of fit test
    virtual std::tuple<bool, GoF> gof(const ModelParams& params, const std::vector<double>& samples) = 0;
};

#endif // ESTIMATOR_H