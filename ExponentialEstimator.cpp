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
#ifdef _OPENMP
    #include <omp.h>
    #pragma message("Compiling ExponentialEstimator with OpenMP support.")
#else
    #pragma message("Compiling ExponentialEstimator without OpenMP support.")
#endif
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

    const size_t len = samples.size();

    // Validate minimum sample size requirement
    if (len < m_min_len) {
        throw invalid_argument("Cannot fit anything with less than 10 values");
    }

    // Calculate initial estimates using method of moments
    const double min = _min(samples);
    const double mean = _mean(samples);
    const double mu = len * (mean - min) / static_cast<double>(len - 1); // Bias-corrected mean estimate

    // Set distribution parameters
    ModelParams params;
    params.alpha = min - mu / static_cast<double>(len);  // Location parameter estimate
    params.beta = 1 / mu;  // Rate parameter (inverse of mean)

    // Perform goodness of fit test
    auto [reject, gof_] = gof(params, samples);

    // Return model only if fit is acceptable
    return reject ? Model() : Model{true, ModelType::EXP, params, gof_};
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

    const size_t len = samples.size();
    const double mu = 1 / params.beta;   // Convert rate to mean for calculations
    const double alpha = params.alpha;   // Location parameter
    const double beta = params.beta;     // Rate parameter
    const double thresh = 1.321;         // Critical value from D'Agostino table 4.14

    // Validate parameters
    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distribution.");
    }
    if (len < m_min_len) {
        throw invalid_argument("Number of samples is not enough or is zero.");
    }

    // Transform data to uniform using exponential CDF
    vector<double> data(len);
    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif
    for(size_t i = 0; i < len; ++i) {
        data[i] = 1.0 - std::exp(-(samples[i] - alpha) / mu);
    }
    sort(data.begin(), data.end());

    // Calculate Anderson-Darling statistic (A²)
    double accum = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(len > 1000)
    #endif    
    for (size_t i = 0; i < len; ++i) {
        const double term1 = (2.0 * (i + 1) - 1.0) * std::log(data[i]);
        const double term2 = (2.0 * len + 1.0 - 2.0 * (i + 1)) * std::log(1.0 - data[i]);
        accum += term1 + term2;
    }

    // Apply correction factors
    const double lenf = static_cast<double>(len);
    const double A2 = (-lenf - accum/lenf) * 
                     (1.0 + 5.4/lenf - 11.0/(lenf * lenf));

    return {A2 > thresh, {A2, thresh}};
}

/**
 * Calculates the cumulative distribution function (CDF) for a single value from 
 * the exponential distribution. Implements the analytical form of the CDF based
 * on the two-parameter exponential distribution.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter (minimum possible value)
 *        - beta: rate parameter (1/mean)
 * @param sample Single input value to evaluate
 * 
 * @return CDF value F(x) = 1 - exp(-beta * (x - alpha))
 * 
 * Implementation details:
 * - Handles boundary conditions implicitly
 * - Uses location-scale parameterization
 * - Based on D'Agostino & Stephens (1986), p. 133
 */
double ExponentialEstimator::cdf(const ModelParams& params, const double& sample) {

    const double alpha = params.alpha;
    const double beta = params.beta;
    const double diff = sample - alpha;
    
    return (diff <= 0) ? 0.0 : 1.0 - exp(-beta * diff);
}

/**
 * Calculates the cumulative distribution function (CDF) for multiple values from
 * the exponential distribution. Provides vectorized implementation for efficiency
 * on large datasets.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter (minimum possible value)
 *        - beta: rate parameter (1/mean)
 * @param samples Vector of input values to evaluate
 * 
 * @return Vector of CDF values F(x) = 1 - exp(-beta * (x - alpha))
 * 
 * Implementation details:
 * - Pre-computes common terms outside the loop for efficiency
 * - Sequential implementation for better cache utilization
 * - Returns vector of same size as input samples
 */
vector<double> ExponentialEstimator::cdf(const ModelParams& params, const vector<double>& samples) {

    const double alpha = params.alpha;
    const double beta = params.beta;

    const size_t len = samples.size();
    vector<double> cdf(len);
    
    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif
    for (size_t i = 0; i < len; ++i) {
        const double diff = samples[i] - alpha;
        cdf[i] = (diff <= 0) ? 0.0 : 1.0 - exp(-beta * diff);
    }
    return cdf;
}

/**
 * Calculates the probability density function (PDF) for a single value from
 * the exponential distribution. Implements the analytical form of the PDF.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter
 *        - beta: rate parameter (1/mean)
 * @param sample Single input value to evaluate
 * 
 * @return PDF value f(x) = beta * exp(-beta * (x - alpha))
 * 
 * Implementation details:
 * - Validates beta parameter (must be positive)
 * - Based on D'Agostino & Stephens (1986), p. 133
 * 
 * @throws invalid_argument if beta <= 0
 */
double ExponentialEstimator::pdf(const ModelParams& params, const double& sample) {

    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distr.");
    }

    const double diff = sample - alpha;
    return (diff <= 0) ? 0.0 : beta * exp(-beta * diff);
}

/**
 * Calculates the probability density function (PDF) for multiple values from
 * the exponential distribution. Provides vectorized implementation for efficiency.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter
 *        - beta: rate parameter (1/mean)
 * @param samples Vector of input values to evaluate
 * 
 * @return Vector of PDF values f(x) = beta * exp(-beta * (x - alpha))
 * 
 * Implementation details:
 * - Validates beta parameter (must be positive)
 * 
 * @throws invalid_argument if beta <= 0
 */
vector<double> ExponentialEstimator::pdf(const ModelParams& params, const vector<double>& samples) {

    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distr.");
    }

    const size_t len = samples.size();
    vector<double> pdf(len);

    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif
    for (size_t i = 0; i < len; ++i) {
        const double diff = samples[i] - alpha;
        pdf[i] = (diff <= 0) ? 0.0 : beta * exp(-beta * diff);
    }
    return pdf;
}

/**
 * Generates a single random value from the exponential distribution using
 * the inverse transform sampling method. Uses the uniform distribution
 * inherited from the base class.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter
 *        - beta: rate parameter (1/mean)
 * @return Random value following the exponential distribution
 * 
 * Uses inverse transform sampling with the formula:
 * X = -ln(1-u)/beta + alpha
 * where u ~ U(0,1)
 * 
 * Implementation details:
 * - Validates beta parameter (must be positive)
 * - Uses uniform distribution from base class
 * - Applies location-scale transformation
 * 
 * @throws invalid_argument if beta <= 0
 */
double ExponentialEstimator::rnd(const ModelParams& params) {

    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distribution");
    }

    const double p = m_unif_dist(m_rnd_gen);
    return -log(1.0 - p) / beta + alpha;
}    

/**
 * Generates multiple random values from the exponential distribution using
 * vectorized inverse transform sampling. Optimized for generating large
 * numbers of random values efficiently.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter
 *        - beta: rate parameter (1/mean)
 * @param length Number of random values to generate
 * 
 * @return Vector of random values following the exponential distribution
 * 
 * Implementation details:
 * - Validates beta parameter (must be positive)
 * 
 * @throws invalid_argument if beta <= 0
 */
vector<double> ExponentialEstimator::rnd(const ModelParams& params, const unsigned& length) {

    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distribution");
    }

    vector<double> samples(length);

    #ifdef _OPENMP
        #pragma omp parallel for if(length > 1000)
    #endif
    for (unsigned i = 0; i < length; ++i) {
        const double p = m_unif_dist(m_rnd_gen);
        samples[i] = -log(1.0 - p) / beta + alpha;
    }
    return samples;
}

/**
 * Calculates the expectation (mean) of the exponential distribution.
 * Implements the analytical formula for the expectation based on the
 * two-parameter exponential distribution.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter (minimum possible value)
 *        - beta: rate parameter (1/mean)
 * @return Expected value E[X] = alpha + 1/beta
 * 
 * Implementation details:
 * - Uses location-scale parameterization
 * - Direct implementation of analytical formula
 */
double ExponentialEstimator::expectation(const ModelParams& params) {
    return params.alpha + 1.0/params.beta;
}

/**
 * Calculates the variance of the exponential distribution.
 * Implements the analytical formula for the variance based on the
 * two-parameter exponential distribution.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter (minimum possible value)
 *        - beta: rate parameter (1/mean)
 * @return Variance Var[X] = 1/beta^2
 * 
 * Implementation details:
 * - Uses rate parameter only (location does not affect variance)
 * - Pre-computes beta^2 for efficiency
 */
double ExponentialEstimator::variance(const ModelParams& params) {
    const double beta = params.beta;
    return 1.0/(beta*beta);
}

/**
 * Calculates the mode of the exponential distribution.
 * The mode is the location parameter since the exponential
 * distribution is monotonically decreasing.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: location parameter (minimum possible value)
 *        - beta: rate parameter (1/mean)
 * @return Mode = alpha (location parameter)
 * 
 * Implementation details:
 * - Mode is independent of rate parameter
 * - Equal to location parameter by definition
 */
double ExponentialEstimator::mode(const ModelParams& params) {
    return params.alpha;
}