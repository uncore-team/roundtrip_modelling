#include "exponential_model.h"
#include <cmath>
#include <stdexcept>

ExponentialModel::ExponentialModel() : lambda(0.0), is_fitted(false) {}

void ExponentialModel::fit(const std::vector<double>& data) {
    if (data.empty()) {
        throw std::invalid_argument("Data cannot be empty.");
    }

    double sum = 0.0;
    for (double value : data) {
        if (value < 0) {
            throw std::invalid_argument("Data values must be non-negative.");
        }
        sum += value;
    }

    lambda = data.size() / sum;
    is_fitted = true;
}

double ExponentialModel::probability(double x) const {
    if (!is_fitted) {
        throw std::runtime_error("Model is not fitted yet.");
    }
    if (x < 0) {
        return 0.0;
    }
    return lambda * std::exp(-lambda * x);
}

double ExponentialModel::get_lambda() const {
    if (!is_fitted) {
        throw std::runtime_error("Model is not fitted yet.");
    }
    return lambda;
}