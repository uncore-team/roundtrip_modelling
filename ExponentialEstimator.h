#ifndef EXPONENTIAL_ESTIMATOR_H
#define EXPONENTIAL_ESTIMATOR_H

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

    /**
     * @brief Calculates the cumulative distribution function (CDF) for a single value.
     * 
     * @param params Distribution parameters {alpha, beta} where:
     *        - alpha: location parameter (minimum possible value)
     *        - beta: rate parameter (1/mean)
     * @param sample Single input value to evaluate
     * @return F(x) = 1 - exp(-beta * (x - alpha))
     */
    virtual double cdf(const ModelParams& params, const double& sample) override;

    /**
     * @brief Calculates the CDF for a vector of values.
     * 
     * @param params Distribution parameters {alpha, beta}
     * @param samples Vector of input values to evaluate
     * @return Vector of CDF values F(x) = 1 - exp(-beta * (x - alpha))
     */
    virtual vector<double> cdf(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Calculates the probability density function (PDF) for a single value.
     * 
     * @param params Distribution parameters {alpha, beta} where:
     *        - alpha: location parameter
     *        - beta: rate parameter (1/mean)
     * @param sample Single input value to evaluate
     * @return f(x) = beta * exp(-beta * (x - alpha))
     * @throws invalid_argument if beta <= 0
     */
    virtual double pdf(const ModelParams& params, const double& sample) override;

    /**
     * @brief Calculates the PDF for a vector of values.
     * 
     * @param params Distribution parameters {alpha, beta}
     * @param samples Vector of input values to evaluate
     * @return Vector of PDF values f(x) = beta * exp(-beta * (x - alpha))
     * @throws invalid_argument if beta <= 0
     */
    virtual vector<double> pdf(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Generates a random value from the exponential distribution.
     * 
     * @param params Distribution parameters {alpha, beta}
     * @return Single random value from the exponential distribution
     * @throws invalid_argument if beta <= 0
     */
    virtual double rnd(const ModelParams& params) override;

    /**
     * @brief Generates multiple random values from the exponential distribution.
     * 
     * @param params Distribution parameters {alpha, beta}
     * @param length Number of random values to generate
     * @return Vector of random values from the exponential distribution
     * @throws invalid_argument if beta <= 0
     */
    virtual vector<double> rnd(const ModelParams& params, const unsigned& length) override;

    /**
     * @brief Calculates the expected value (mean) of the exponential distribution
     * 
     * @param params Distribution parameters {alpha, beta} where:
     *        - alpha: location parameter
     *        - beta: rate parameter (1/mean)
     * @return Expected value E[X] = alpha + 1/beta
     */
    virtual double expectation(const ModelParams& params) override;

    /**
     * @brief Calculates the variance of the exponential distribution
     * 
     * @param params Distribution parameters {alpha, beta} where:
     *        - alpha: location parameter
     *        - beta: rate parameter (1/mean)
     * @return Variance Var[X] = 1/beta^2
     */
    virtual double variance(const ModelParams& params) override;

    /**
     * @brief Calculates the mode of the exponential distribution
     * 
     * @param params Distribution parameters {alpha, beta} where:
     *        - alpha: location parameter
     *        - beta: rate parameter (1/mean)
     * @return Mode = alpha (location parameter)
     */
    virtual double mode(const ModelParams& params) override;
};

#endif // EXPONENTIAL_ESTIMATOR_H