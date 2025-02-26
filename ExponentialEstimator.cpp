/**
 * @brief Implementation of the ExponentialEstimator class for two-parameter exponential distribution.
 * 
 * Provides methods for:
 * - Parameter estimation using maximum likelihood
 * - Goodness of fit testing using Anderson-Darling
 * - Distribution fitting with validation
 * 
 * Based on D'Agostino & Stephens (1986), Chapter 4.
 */

#include "ExponentialEstimator.h"

using namespace std;

/**
 * Constructor for ExponentialEstimator class.
 * Initializes the base Estimator class with a minimum required sample size of 10.
 * This minimum size ensures reliable parameter estimation for the exponential distribution.
 */
ExponentialEstimator::ExponentialEstimator() : Estimator(10) {
}

/**
 * Fits a two-parameter exponential distribution to the given samples.
 * Implementation based on D'Agostino & Stephens (1986), p. 141.
 * 
 * @param samples Input vector of observed values
 * 
 * @return Model structure containing:
 *         - defined: true if fit was successful
 *         - type: ModelType::EXP for successful fit, ModelType::None otherwise
 *         - params: {alpha, beta} parameters where:
 *           * alpha: location parameter (minimum possible value)
 *           * beta: rate parameter (1/mean)
 *         - gof: Goodness of fit statistics
 * 
 * Implementation details:
 * - Corrected implementation from D'Agostino's book
 * - beta is calculated as 1/mean (reciprocal of the mean)
 * - Uses maximum likelihood estimation for parameters
 * - Performs Anderson-Darling goodness of fit test
 * 
 * @throws invalid_argument if samples.size() < m_min_len
 */
Model ExponentialEstimator::fit(const vector<double>& samples) {
    int len = samples.size();

    // Validate minimum sample size requirement
    if (len < m_min_len) {
        throw invalid_argument("Cannot fit anything with less than 10 values");
    }

    // Calculate initial estimates using method of moments
    double min = _min(samples);
    double mean = _mean(samples);
    double mu = len * (mean - min) / (len - 1); // Bias-corrected mean estimate

    // Set distribution parameters
    ModelParams params;
    params.alpha = min - mu / len;    // Location parameter estimate
    params.beta = 1 / mu;             // Rate parameter (inverse of mean)

    // Perform goodness of fit test
    auto [reject, gof_] = gof(params, samples);

    // Return model only if fit is acceptable
    if (!reject) {
        return {true, ModelType::EXP, params, gof_};
    }
    else { 
        return Model(); // Return invalid model if fit is rejected
    }
}

/**
 * Performs Anderson-Darling goodness of fit test for exponential distribution.
 * Based on D'Agostino & Stephens (1986), p. 141, Case 3 (both parameters unknown).
 * 
 * @param params Model parameters {alpha, beta} to test
 * @param samples Input vector of observed values
 * 
 * @return tuple<bool, GoF> containing:
 *         - bool: true if null hypothesis should be rejected (poor fit)
 *         - GoF: {statistic, threshold} where:
 *           * statistic: modified A² test statistic
 *           * threshold: critical value (1.321) from Table 4.14
 * 
 * Implementation details:
 * - Transforms data to uniform distribution using CDF
 * - Applies small sample correction factor: (1 + 5.4/n - 11/n²)
 * - Uses right-tailed test with significance level 0.05
 * - Rejection criteria: A² > 1.321
 * 
 * @throws invalid_argument if:
 *         - beta <= 0
 *         - samples.size() < m_min_len
 */
tuple<bool, GoF> ExponentialEstimator::gof(const ModelParams& params, const vector<double>& samples) {
    int len = samples.size();
    double mu = 1 / params.beta;      // Convert rate to mean for calculations
    double alpha = params.alpha;       // Location parameter
    double beta = params.beta;         // Rate parameter
    double thresh = 1.321;            // Critical value from D'Agostino table 4.14

    // Validate parameters
    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distribution.");
    }
    if (len < m_min_len) {
        throw invalid_argument("Number of samples is not enough or is zero.");
    }

    // Transform data to uniform using exponential CDF
    vector<double> data(len);
    transform(
        samples.begin(), samples.end(),
        data.begin(),
        [alpha, mu](double sample) { return 1 - exp(-(sample - alpha) / mu); }
    );
    sort(data.begin(), data.end());

    // Calculate Anderson-Darling statistic (A²)
    double accum = 0.0;
    for (int i = 0; i < len; ++i) {
        accum += (2 * (i + 1) - 1) * log(data[i]) + 
                 (2 * len + 1 - 2 * (i + 1)) * log(1 - data[i]);
    }

    // Compute test statistic with sample size correction
    double A2 = -len - (1.0 / len) * accum;
    A2 *= (1 + 5.4 / len - 11.0 / (len * len));  // Small sample correction

    // Compare against critical value
    bool reject = (A2 > thresh);  // Reject if statistic exceeds threshold

    return {reject, {A2, thresh}};
}
