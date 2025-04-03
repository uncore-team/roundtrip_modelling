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
#include <interpolation.h>
#ifdef _OPENMP
    #include <omp.h>
    //#pragma message("Compiling ExponentialEstimator with OpenMP support.")
#else
    #pragma message("Compiling ExponentialEstimator without OpenMP support.")
#endif
#include "ExponentialEstimator.h"

using namespace std;
using namespace alglib;

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

    // PROBLEM: This method produces negative alpha in a high proportion of
    // the estimations
    // const double min = _min(samples);
    // const double mean = _mean(samples);
    // const double mu = len*(mean - min)/static_cast<double>(len - 1); // Bias-corrected mean estimate

    // // Set distribution parameters
    // ModelParams params;
    // params.alpha = min - mu/static_cast<double>(len);  // Location parameter estimate
    // params.beta = 1/mu;  // Rate parameter (inverse of mean)

    // Calculate parameters using MLE estimation that forces 
    // alpha to be positive and beta to be finite and positive
    const double min = _min(samples);

    double accum = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(len > OMP_THRESH)
    #endif    
    for (const double sample: samples) {
        accum += (sample - min);
    }
    const double mean = accum/static_cast<double>(len);

    // Set distribution parameters
    ModelParams params;
    params.alpha = min;  // Location parameter estimate
    params.beta = 1/mean;  // Rate parameter (inverse of mean)

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
tuple<bool, GoF> ExponentialEstimator::gof(const ModelParams& params, const vector<double>& samples, bool prev_model) {

    // Parameters
    const double alpha = params.alpha;   // Location parameter
    const double beta = params.beta;     // Rate parameter

    // Validate parameters
    if (alpha < 0.0) {
        throw invalid_argument("Invalid Exponential distr.: alpha < 0");
    }
    if (beta <= 0.0) {
        throw invalid_argument("Invalid Exponential distr.: beta <= 0");
    }

    const double min = _min(samples);
    if (min < alpha) {
        return {true, {Inf, NaN}}; // cannot accept a distribution if some value falls below its minimum
    }

    // Transform data to uniform using exponential CDF
    const size_t len = samples.size();
    vector<double> data(len);
    #ifdef _OPENMP
        #pragma omp parallel for if(len > OMP_THRESH)
    #endif
    for (size_t i = 0; i < len; ++i) {
        data[i] = 1.0 - exp(-beta*(samples[i] - alpha));
    }    
    sort(data.begin(), data.end());

    // Calculate Anderson-Darling statistic (A²) This statistic measures the 
    // squared distance between experimental and theoretical Zs, and, indirectly,
    // between theoretical and experimental Xs (p. 100)
    const double lenf = static_cast<double>(len);
    double accum = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(len > OMP_THRESH)
    #endif
    for (size_t i = 0; i < len; ++i) {
        accum += (2.0*(i + 1) - 1.0)*log(data[i]) + (2.0*lenf + 1.0 - 2.0*(i + 1))*log(1.0 - data[i]);
    }
    double A2 = -lenf - accum/lenf;

    // Apply correction factors
    double thresh;
    if (prev_model) {
    	thresh = 2.492; // threshold for the case that parameters do not come from sample (n >= 5)
                        // D'Agostino & Stephens (1986), 0.05 significance level, p. 105, table 4.2, right tail
    }
    else {
        if (len > 100) {
            thresh = 1.321;  // threshold for the case of parameters coming from sample (case 3)
                             // D'Agostino & Stephens (1986), 0.05 significance level, p. 141, table 4.14, right tail
        }
        else {
            real_1d_array ns = "[5, 10, 15, 20, 25, 50, 100]"; // lengths listed in D'Agostino table 4.15
            real_1d_array ts = "[0.725, 0.920, 1.009, 1.062, 1.097, 1.197, 1.250]";
    
            spline1dinterpolant s;
            spline1dbuildcatmullrom(ns, ts, s); // best interpolation since points are non-linear
                                                // do not include the point at > 100 since it is
                                                // too far and distorts the spline strongly
            thresh = spline1dcalc(s, lenf);
        }
        A2 *= (1.0 + 5.4/lenf - 11.0/(lenf*lenf)); // apply correction factors
    }
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

    // Parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0.0) {
        throw invalid_argument("Invalid beta for exponential distr.");
    }

    const double z = sample - alpha;
    return (z <= 0.0) ? 0.0 : 1.0 - exp(-beta*z);
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

    // Parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distr.");
    }

    // Calculate the CDF of the exponential distribution
    const size_t len = samples.size();
    vector<double> cdf(len);

    #ifdef _OPENMP
        #pragma omp parallel for if(len > OMP_THRESH)
    #endif
    for (size_t i = 0; i < len; ++i) {
        const double z = samples[i] - alpha;
        cdf[i] = (z <= 0.0) ? 0.0 : 1.0 - exp(-beta*z);
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

    // Parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distr.");
    }

    const double z = sample - alpha;
    return (z < 0.0) ? 0.0 : beta*exp(-beta*z);
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

    // Parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distr.");
    }

    // Calculate the PDF of the exponential distribution
    const size_t len = samples.size();
    vector<double> pdf(len);

    #ifdef _OPENMP
        #pragma omp parallel for if(len > OMP_THRESH)
    #endif
    for (size_t i = 0; i < len; ++i) {
        const double z = samples[i] - alpha;
        pdf[i] = (z < 0.0) ? 0.0 : beta*exp(-beta*z);
    }
    return pdf;
}

/**
 * Generates a single random value from the exponential distribution using
 * the C++ Standard Library's exponential_distribution.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: Location parameter (shift)
 *        - beta: Rate parameter (1/mean, must be positive)
 * @return Random value following the exponential distribution
 * 
 * Implementation details:
 * - Uses exponential_distribution for exponential random number generation
 * - Validates beta parameter (must be positive)
 * - Applies location-scale transformation: X = alpha + Y, where Y ~ Exp(beta)
 * 
 * @throws invalid_argument if beta <= 0
 */
double ExponentialEstimator::rnd(const ModelParams& params) {

    // Extract parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    // Validate beta parameter
    if (beta <= 0.0) {
        throw invalid_argument("Invalid beta for exponential distribution (must be > 0)");
    }

    // Create an exponential distribution with rate parameter beta
    exponential_distribution<double> exponential(beta);

    // Generate a random value from the exponential distribution and apply the shift alpha
    return alpha + exponential(m_rnd_gen);
}

/**
 * Generates multiple random values from the exponential distribution using
 * the C++ Standard Library's exponential_distribution.
 * 
 * @param params Distribution parameters {alpha, beta} where:
 *        - alpha: Location parameter (shift)
 *        - beta: Rate parameter (1/mean, must be positive)
 * @param length Number of random values to generate
 * @return Vector of random values following the exponential distribution
 * 
 * Implementation details:
 * - Uses exponential_distribution for exponential random number generation
 * - Validates beta parameter (must be positive)
 * - Applies location-scale transformation: X = alpha + Y, where Y ~ Exp(beta)
 * - Uses OpenMP for parallel generation when length > OMP_THRESH
 * 
 * @throws invalid_argument if beta <= 0
 */
vector<double> ExponentialEstimator::rnd(const ModelParams& params, const unsigned& length) {

    // Extract parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    // Validate beta parameter
    if (beta <= 0.0) {
        throw invalid_argument("Invalid beta for exponential distribution (must be > 0)");
    }

    // Create an exponential distribution with rate parameter beta
    exponential_distribution<double> exponential(beta);

    // Generate 'length' random values from the exponential distribution
    vector<double> rnd(length);

    #ifdef _OPENMP
        #pragma omp parallel for if(length > OMP_THRESH)
    #endif
    for (size_t i = 0; i < length; ++i) {
        rnd[i] = alpha + exponential(m_rnd_gen);
    }

    return rnd;
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

   // Parameters
    const double alpha = params.alpha;
    const double beta = params.beta;

    return alpha + 1.0/beta;
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

    // Parameters
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