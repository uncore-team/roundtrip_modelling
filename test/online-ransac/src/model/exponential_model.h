#ifndef EXPONENTIAL_MODEL_H
#define EXPONENTIAL_MODEL_H

#include <vector>

class ExponentialModel {
public:
    ExponentialModel();

    void fit(const std::vector<double>& data);
    double predict(double x) const;
    double getLambda() const;

private:
    double lambda; // Rate parameter for the exponential distribution
    bool isFitted; // Flag to check if the model has been fitted
};

#endif // EXPONENTIAL_MODEL_H