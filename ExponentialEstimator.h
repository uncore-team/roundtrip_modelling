#ifndef EXPONENTIAL_ESTIMATOR_H
#define EXPONENTIAL_ESTIMATOR_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Estimator.h"

using namespace std;

/**
 * @brief Class for fitting and analyzing two-parameter exponential distributions.
 * 
 * Implements maximum likelihood estimation and goodness-of-fit testing for
 * the exponential distribution with parameters:
 * - alpha (location)
 * - beta (rate)
 * 
 * The probability density function is:
 * f(x) = beta * exp(-beta*(x - alpha)) for x ≥ alpha
 * 
 * Implementation based on D'Agostino & Stephens (1986), Chapter 4.
 */
class ExponentialEstimator : public Estimator {
public:
    /**
     * @brief Constructor initializing minimum sample size requirement.
     * Sets minimum sample size to 10 for reliable parameter estimation.
     */
    ExponentialEstimator();

    /**
     * @brief Fits a two-parameter exponential distribution to sample data.
     * 
     * @param samples Vector of observations to fit
     * @return Model structure containing:
     *         - defined: true if fit succeeded
     *         - type: ModelType::EXP or ModelType::None
     *         - params: {alpha, beta} distribution parameters
     *         - gof: Goodness of fit statistics
     * @throws invalid_argument if samples.size() < 10
     */
    Model fit(const vector<double>& samples) override;

    /**
     * @brief Performs Anderson-Darling goodness-of-fit test.
     * 
     * @param params Distribution parameters {alpha, beta} to test
     * @param samples Vector of observations to test against
     * @return tuple<bool, GoF>:
     *         - bool: true if fit should be rejected
     *         - GoF: {test statistic, critical value}
     * @throws invalid_argument if beta <= 0 or samples.size() < 10
     */
    tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) override;
};

#endif // EXPONENTIAL_ESTIMATOR_H