#ifndef LOGLOGISTIC_ESTIMATOR_H
#define LOGLOGISTIC_ESTIMATOR_H

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
    tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples, bool prev_model = false) override;  // bool previous_model

    /**
     * @brief Calculates the CDF for a single value from the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c} where:
     *        - a: location parameter (minimum possible value)
     *        - b: scale parameter
     *        - c: shape parameter
     * @param sample Single input value to evaluate
     * @return CDF value F(x) = 1/(1 + ((x-a)/b)^(-1/c))
     */
    virtual double cdf(const ModelParams& params, const double& sample) override;

    /**
     * @brief Calculates the CDF for multiple values from the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c}
     * @param samples Vector of input values to evaluate
     * @return Vector of CDF values F(x) = 1/(1 + ((x-a)/b)^(-1/c))
     */
    virtual vector<double> cdf(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Calculates the PDF for a single value from the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c}
     * @param sample Single input value to evaluate
     * @return PDF value f(x) = (c/b)((x-a)/b)^(c-1)/(1+((x-a)/b)^c)^2
     */
    virtual double pdf(const ModelParams& params, const double& sample) override;

    /**
     * @brief Calculates the PDF for multiple values from the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c}
     * @param samples Vector of input values to evaluate
     * @return Vector of PDF values
     */
    virtual vector<double> pdf(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Generates a single random value from the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c}
     * @return Random value following the log-logistic distribution
     */
    virtual double rnd(const ModelParams& params) override;

    /**
     * @brief Generates multiple random values from the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c}
     * @param length Number of random values to generate
     * @return Vector of random values following the log-logistic distribution
     */
    virtual vector<double> rnd(const ModelParams& params, const unsigned& length) override;

    /**
     * @brief Calculates the expected value (mean) of the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c} where:
     *        - a: location parameter
     *        - b: scale parameter
     *        - c: shape parameter
     * @return Expected value E[X] when it exists (c > 1)
     */
    virtual double expectation(const ModelParams& params) override;

    /**
     * @brief Calculates the variance of the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c}
     * @return Variance Var[X] when it exists (c > 2)
     */
    virtual double variance(const ModelParams& params) override;

    /**
     * @brief Calculates the mode (most frequent value) of the log-logistic distribution
     *
     * @param params Distribution parameters {a, b, c} where:
     *        - a: location parameter
     *        - b: scale parameter
     *        - c: shape parameter
     * @return Mode = a + b*((c-1)/(c+1))^(1/c) when c > 1, a otherwise
     */
    virtual double mode(const ModelParams& params) override;

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
