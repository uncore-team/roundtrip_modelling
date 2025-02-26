#ifndef LOGLOGISTIC_ESTIMATOR_H
#define LOGLOGISTIC_ESTIMATOR_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Estimator.h"

using namespace std;

/**
 * @brief Class for fitting and analyzing three-parameter log-logistic distributions.
 * 
 * Implements maximum likelihood estimation and goodness-of-fit testing for
 * the log-logistic distribution with parameters:
 * - a (location)
 * - b (scale)
 * - c (shape)
 * 
 * The probability density function is:
 * f(x) = (c/b)*((x-a)/b)^(c-1)/(1+((x-a)/b)^c)² for x > a
 * 
 * Implementation uses a two-stage optimization approach:
 * 1. Shape parameter estimation with fixed location and scale
 * 2. Full three-parameter optimization using Levenberg-Marquardt
 */
class LogLogisticEstimator : public Estimator {
public:
    /**
     * @brief Constructor initializing minimum sample size requirement.
     * Sets minimum sample size to 10 for reliable parameter estimation.
     */
    LogLogisticEstimator();

    /**
     * @brief Fits a three-parameter log-logistic distribution to sample data.
     * 
     * @param samples Vector of observations to fit
     * @return Model structure containing:
     *         - defined: true if fit succeeded
     *         - type: ModelType::LL3 or ModelType::None
     *         - params: {a, b, c} distribution parameters
     *         - gof: Goodness of fit statistics
     * @throws invalid_argument if samples.size() < 10
     */
    Model fit(const vector<double>& samples) override;

    /**
     * @brief Performs Anderson-Darling goodness-of-fit test.
     * 
     * @param params Distribution parameters {a, b, c} to test
     * @param samples Vector of observations to test against
     * @return tuple<bool, GoF>:
     *         - bool: true if fit should be rejected
     *         - GoF: {test statistic, critical value}
     * @throws invalid_argument if b <= 0 or c <= 0 or samples.size() < 10
     */
    tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) override;  // bool previous_model

private:
    /**
     * @brief Estimates shape parameter c keeping location and scale fixed.
     * 
     * @param samples Input data vector
     * @param a Fixed location parameter
     * @param b Fixed scale parameter
     * @param c Input/Output shape parameter
     * @return Optimization termination type (>0 success, ≤0 failure)
     */
    int c_estimate( const vector<double>& samples, const double& a, const double& b, double& c);

    /**
     * @brief Estimates all three parameters simultaneously.
     * 
     * @param samples Input data vector
     * @param a Input/Output location parameter
     * @param b Input/Output scale parameter
     * @param c Input/Output shape parameter
     * @return Optimization termination type (>0 success, ≤0 failure)
     */    
    int abc_estimate(const vector<double>& samples, double& a, double& b, double& c);
};

#endif // LOGLOGISTIC_ESTIMATOR_H