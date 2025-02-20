// filepath: /ransaconline2/ransaconline2/src/estimators/LL3Estimator.h
#ifndef LL3ESTIMATOR_H
#define LL3ESTIMATOR_H

#include <vector>
#include <optional>
#include "Model.h"

// LL3Estimator class for fitting a log-logistic model to data.
class LL3Estimator {
public:
    // Constructor
    LL3Estimator();

    // Fit the log-logistic model to the provided data.
    // Parameters:
    // - data: A vector of double values representing the data to fit.
    // Returns: A boolean indicating success or failure of the fitting process.
    bool fit(const std::vector<double>& data);

    // Assess the goodness of fit of the model to the data.
    // Parameters:
    // - data: A vector of double values representing the data to assess.
    // Returns: A tuple containing a boolean indicating rejection of the null hypothesis,
    //          the statistic of the goodness of fit test, and the threshold value.
    std::tuple<bool, double, double> assess(const std::vector<double>& data) const;

    // Get the fitted model.
    // Returns: An optional Model object containing the fitted model.
    std::optional<Model> getModel() const;

private:
    Model model; // The log-logistic model fitted to the data.
};

#endif // LL3ESTIMATOR_H