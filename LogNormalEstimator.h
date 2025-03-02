#ifndef LOGNORMAL_ESTIMATOR_H
#define LOGNORMAL_ESTIMATOR_H

#include "Estimator.h"

using namespace std;

// LogNormalEstimator class for fitting and assessing an LogLogistic3 model to data.
class LogNormalEstimator : public Estimator {
public:
    /**
     * @brief Constructor initializing minimum sample size requirement.
     * Sets minimum sample size to 10 for reliable parameter estimation.
     */
    LogNormalEstimator();

    /**
     * @brief Fits the log-normal distribution parameters to the provided sample data
     * using maximum likelihood estimation.
     * 
     * @param samples Vector of observations to fit the distribution to
     * @return Model containing the fitted parameters {μ, σ} where:
     *         - μ: location parameter
     *         - σ: scale parameter (σ > 0)
     * 
     * Implementation details:
     * - Uses method of moments for initial parameter estimates
     * - Optimizes likelihood using numerical methods
     * - Validates parameter constraints (σ > 0)
     */
    Model fit(const vector<double>& samples) override;

    /**
     * @brief Performs Anderson-Darling goodness of fit test for the log-normal distribution
     * with the given parameters.
     * 
     * @param params Distribution parameters {μ, σ} where:
     *        - μ: location parameter
     *        - σ: scale parameter (σ > 0)
     * @param samples Vector of observations to test against the distribution
     * @return Tuple containing:
     *         - bool: true if null hypothesis is rejected
     *         - GoF: {statistic, critical value} for the test
     * 
     * Implementation details:
     * - Uses Anderson-Darling test statistic
     * - Compares against critical values for log-normal distribution
     * - Requires minimum sample size for validity
     */
    tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Calculates the cumulative distribution function (CDF) for a single value
     * from the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ} where:
     *        - μ: location parameter
     *        - σ: scale parameter (σ > 0)
     * @param sample Single input value to evaluate
     * @return CDF value F(x) = Φ((ln(x) - μ)/σ) where Φ is the standard normal CDF
     */
    virtual double cdf(const ModelParams& params, const double& sample) override;

    /**
     * @brief Calculates the cumulative distribution function (CDF) for multiple values
     * from the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @param samples Vector of input values to evaluate
     * @return Vector of CDF values F(x) = Φ((ln(x) - μ)/σ)
     */
    virtual vector<double> cdf(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Calculates the probability density function (PDF) for a single value
     * from the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @param sample Single input value to evaluate
     * @return PDF value f(x) = 1/(xσ√(2π)) * exp(-(ln(x)-μ)²/(2σ²))
     */
    virtual double pdf(const ModelParams& params, const double& sample) override;

    /**
     * @brief Calculates the probability density function (PDF) for multiple values
     * from the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @param samples Vector of input values to evaluate
     * @return Vector of PDF values
     */
    virtual vector<double> pdf(const ModelParams& params, const vector<double>& samples) override;

    /**
     * @brief Generates a single random value from the log-normal distribution using
     * the inverse transform sampling method.
     * 
     * @param params Distribution parameters {μ, σ}
     * @return Random value X = exp(μ + σZ) where Z ~ N(0,1)
     */
    virtual double rnd(const ModelParams& params) override;

    /**
     * @brief Generates multiple random values from the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @param length Number of random values to generate
     * @return Vector of random values following the log-normal distribution
     */
    virtual vector<double> rnd(const ModelParams& params, const unsigned& length) override;

    /**
     * @brief Calculates the expected value (mean) of the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @return Expected value E[X] = exp(μ + σ²/2)
     */
    virtual double expectation(const ModelParams& params) override;

    /**
     * @brief Calculates the variance of the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @return Variance Var[X] = [exp(σ²) - 1]exp(2μ + σ²)
     */
    virtual double variance(const ModelParams& params) override;

    /**
     * @brief Calculates the mode of the log-normal distribution.
     * 
     * @param params Distribution parameters {μ, σ}
     * @return Mode = exp(μ - σ²)
     */
    virtual double mode(const ModelParams& params) override;
};

#endif // LOGNORMAL_ESTIMATOR_H