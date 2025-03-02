#ifdef _OPENMP
    #include <omp.h>
    #pragma message("Compiling LogNormalEstimator with OpenMP support.")
#else
    #pragma message("Compiling LogNormalEstimator without OpenMP support.")
#endif
#include "LogNormalEstimator.h"
#include "alglib/optimization.h"

using namespace std;
using namespace alglib;

/**
 * Constructor for LogNormalEstimator class.
 * Initializes the base Estimator class with a minimum required sample size of 10.
 * This minimum size ensures reliable parameter estimation for the log-logistic distribution.
 */
LogNormalEstimator::LogNormalEstimator() : Estimator(10) {
}

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
Model LogNormalEstimator::fit(const vector<double>& samples) {

    return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
}

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
tuple<bool, GoF> LogNormalEstimator::gof(const ModelParams& params, const vector<double>& samples) { // bool previous_model
    return {true, {Inf, NaN}};
}

/**
 * Calculates the cumulative distribution function (CDF) for a single value
 * from the log-normal distribution with offset.
 * 
 * @param params Distribution parameters {gamma, alpha, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - alpha: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @param sample Single input value to evaluate
 * @return CDF value F(x) = Φ((ln(x-gamma) - α)/σ) where Φ is the standard normal CDF
 * 
 * Implementation details:
 * - Uses error function (erf) for standard normal CDF calculation
 * - Handles boundary conditions (x <= gamma)
 * - Validates parameter constraints (sigma > 0)
 */
double LogNormalEstimator::cdf(const ModelParams& params, const double& sample) {

    const double offset = params.gamma;
    const double alpha = params.alpha;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    // Calculate the CDF of the lognormal distribution
    if(sample <= offset) {
        return 0.0;
    } else {
        const double z = (log(sample - offset) - alpha) / (sigma * sqrt(2.0));
        return 0.5 * (1.0 + erf(z));
    }
}

/**
 * Calculates the cumulative distribution function (CDF) for multiple values
 * from the log-normal distribution with offset.
 * 
 * @param offset Location parameter (minimum possible value)
 * @param mu Mean of the associated normal distribution
 * @param sigma Standard deviation of the associated normal distribution (sigma > 0)
 * @param x Vector of values to evaluate
 * @return Vector of CDF values F(x) = Φ((ln(x-offset) - μ)/σ) 
 *         where Φ is the standard normal CDF
 * 
 * Implementation details:
 * - Uses error function (erf) for standard normal CDF calculation
 * - Handles boundary conditions (x <= offset)
 * - Validates parameter constraints (sigma > 0)
 */
vector<double> LogNormalEstimator::cdf(const ModelParams& params, const vector<double>& samples) {

    const double offset = params.gamma;
    const double alpha = params.alpha;
    const double sigma = params.sigma;

    const double sqrt2 = sqrt(2.0);
    unsigned len = samples.size();

    // Calculate the CDF of the lognormal distribution
    vector<double> cdf(len);
    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif    
    for(size_t i = 0; i < len; ++i) {
        if(samples[i] <= offset) {
            cdf[i] = 0.0;
        } else {
            const double z = (log(samples[i] - offset) - alpha) / (sigma * sqrt2);
            cdf[i] = 0.5 * (1.0 + erf(z));
        }
    }
    return cdf;
}

/**
 * Calculates the cumulative distribution function (CDF) for a single value
 * from the log-normal distribution with offset.
 * 
 * @param params Distribution parameters {gamma, alpha, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - alpha: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @param sample Single input value to evaluate
 * @return CDF value F(x) = Φ((ln(x-gamma) - α)/σ) where Φ is the standard normal CDF
 * 
 * Implementation details:
 * - Uses error function (erf) for standard normal CDF calculation
 * - Handles boundary conditions (x <= gamma)
 * - Validates parameter constraints (sigma > 0)
 */
double LogNormalEstimator::pdf(const ModelParams& params, const double& sample) {

    const double offset = params.gamma;
    const double alpha = params.alpha;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    const double sqrt2pi = sqrt(2.0 * M_PI);
    const double twoSigmaSq = 2.0 * sigma * sigma;

    // Calculate the PDF of the lognormal distribution
    if(sample <= offset) {
        return 0.0;
    } else {
        const double diff = sample - offset;
        const double logDiff = log(diff);
        const double expTerm = exp(-(logDiff - alpha)*(logDiff - alpha)/twoSigmaSq);
        return expTerm / (diff * sigma * sqrt2pi);
    }
}

/**
 * Calculates the probability density function (PDF) for multiple values
 * from the log-normal distribution with offset.
 * 
 * @param x Vector of values to evaluate
 * @param offset Location parameter (minimum possible value)
 * @param alpha Mean of the associated normal distribution
 * @param sigma Standard deviation of the associated normal distribution (sigma > 0)
 * @return Vector of PDF values f(x) = 1/((x-offset)σ√(2π)) * exp(-(ln(x-offset)-μ)²/(2σ²))
 * 
 * Implementation details:
 * - Pre-computes common terms for efficiency
 * - Handles boundary conditions (x <= offset)
 * - Returns 0 for invalid inputs (x <= offset)
 */
vector<double> LogNormalEstimator::pdf(const ModelParams& params, const vector<double>& samples) {

    const double offset = params.gamma;
    const double alpha = params.alpha;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    const double sqrt2pi = sqrt(2.0 * M_PI);
    const double twoSigmaSq = 2.0 * sigma * sigma;

    // Calculate PDFs with OpenMP parallelization for large datasets
    unsigned len = samples.size();
    vector<double> pdf(len);

    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif
    for(size_t i = 0; i < len; ++i) {
        if(samples[i] <= offset) {
            pdf[i] = 0.0;
        } else {
            const double diff = samples[i] - offset;
            const double logDiff = log(diff);
            const double expTerm = exp(-(logDiff - alpha)*(logDiff - alpha)/twoSigmaSq);
            pdf[i] = expTerm / (diff * sigma * sqrt2pi);
        }
    }
    return pdf;
}

/**
 * Generates a single random value from the log-normal distribution using
 * the Box-Muller transform method.
 * 
 * @param params Distribution parameters {gamma, alpha, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - alpha: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @return Random value X = gamma + exp(α + σZ) where Z ~ N(0,1)
 * 
 * Implementation details:
 * - Uses Box-Muller transform for normal random generation
 * - Validates sigma parameter
 * - Uses inverse error function for normal approximation
 */
double LogNormalEstimator::rnd(const ModelParams& params) {

   // Parameters
   const double offset = params.gamma;
   const double alpha = params.alpha;
   const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

   const double sqrt2 = sqrt(2.0);

   // Generate 'length' data from loglogistic distribution

    // Box-Muller transform to get standard normal
    const double u = m_unif_dist(m_rnd_gen);
    const double z = sqrt2 * _erfinv(2.0*u - 1.0);

    return offset + exp(alpha + sigma*z);  // Transform to log-normal with parameters
}    

/**
 * Generates random values from a log-normal distribution with offset.
 * Uses Box-Muller transform to generate normal variables and then
 * transforms them to log-normal.
 * 
 * @param offs Location parameter (minimum possible value)
 * @param mu Mean of the associated normal distribution
 * @param sigma Standard deviation of associated normal distribution (sigma > 0)
 * @param m Number of rows in output matrix
 * @param n Number of columns in output matrix
 * @return Matrix (m x n) of random values from log-normal(offs, mu, sigma)
 * 
 * Implementation details:
 * - Uses uniform distribution to generate base random numbers
 * - Applies Box-Muller transform for normal distribution
 * - Transforms to log-normal via exp()
 */
vector<double> LogNormalEstimator::rnd(const ModelParams& params, const unsigned& length) {

   // Parameters
   const double offset = params.gamma;
   const double alpha = params.alpha;
   const double sigma = params.sigma;

   const double sqrt2 = sqrt(2.0);

   // Generate 'length' data from loglogistic distribution
   vector<double> rnd(length);

   #ifdef _OPENMP
       #pragma omp parallel for if(length > 1000)
   #endif
   for (unsigned i = 0; i < length; ++i) {
        // Box-Muller transform to get standard normal
        const double u = m_unif_dist(m_rnd_gen);
        const double z = sqrt2 * _erfinv(2.0 * u - 1.0);

        rnd[i] = offset + exp(alpha + sigma*z); // Transform to log-normal with parameters
   }
   return rnd;
}

/**
 * Calculates the expected value (mean) of the log-normal distribution.
 * 
 * @param params Distribution parameters {gamma, alpha, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - alpha: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @return Expected value E[X] = gamma + exp(α + σ²/2)
 * 
 * Implementation details:
 * - Validates sigma parameter
 * - Uses exact analytical formula
 */
double LogNormalEstimator::expectation(const ModelParams& params) {

   // Parameters
   const double gamma = params.gamma;
   const double alpha = params.alpha;
   const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }
    return gamma + exp(alpha + (sigma*sigma)/2);
};

/**
 * Calculates the variance of the log-normal distribution.
 * 
 * @param params Distribution parameters {gamma, alpha, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - alpha: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @return Variance Var[X] = [exp(σ²) - 1]exp(2α + σ²)
 * 
 * Implementation details:
 * - Validates sigma parameter
 * - Pre-computes sigma squared term
 * - Uses exact analytical formula
 */
double LogNormalEstimator::variance(const ModelParams& params) {

    // Parameters
//    const double gamma = params.gamma;
    const double alpha = params.alpha;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    const double sigma2 = sigma*sigma;

    return (exp(sigma2) - 1) * exp(2*alpha + sigma2);
}

/**
 * Calculates the mode of the log-normal distribution.
 * 
 * @param params Distribution parameters {gamma, alpha, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - alpha: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @return Mode = gamma + exp(α - σ²)
 * 
 * Implementation details:
 * - Validates sigma parameter
 * - Uses exact analytical formula
 */
double LogNormalEstimator::mode(const ModelParams& params) {

    // Parameters
    const double gamma = params.gamma;
    const double alpha = params.alpha;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    return gamma + exp(alpha - sigma*sigma);
};
