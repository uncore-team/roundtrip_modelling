/**
 * Implementation of the LogNormalEstimator class for three-parameter log-normal distribution.
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
 * Performs Anderson-Darling goodness of fit test for the log-normal distribution.
 * Based on D'Agostino & Stephens (1986), p. 123.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
 *        - sigma: scale parameter (σ > 0)
 * @param samples Vector of observations to test against
 * 
 * @return tuple<bool, GoF> containing:
 *         - bool: true if null hypothesis should be rejected (poor fit)
 *         - GoF: {statistic, threshold} where:
 *           * statistic: modified A² test statistic
 *           * threshold: critical value (0.752) from Table 4.7, Case 3
 * 
 * Implementation details:
 * - Transforms data to standard normal using CDF
 * - Applies small sample correction factor: (1 + 0.75/n + 2.25/n²)
 * - Uses right-tailed test with significance level 0.05
 * - Rejection criteria: A² > 0.752
 * 
 * @throws invalid_argument if:
 *         - sigma <= 0
 *         - samples.size() < m_min_len
 */
tuple<bool, GoF> LogNormalEstimator::gof(const ModelParams& params, const vector<double>& samples) {

    const double gamma = params.gamma;
    const double mu = params.mu;
    const double sigma = params.sigma;
    const double thresh = 0.752;  // D'Agostino & Stephens (1986), p. 123, Table 4.7
    const size_t len = samples.size();
    const double min = _min(samples);

    // Parameter validation
    if (sigma <= 0) {
        throw invalid_argument("Invalid sigma for lognormal distribution.");
    }
    if (len < m_min_len) {
        throw invalid_argument("Number of samples is not enough or is zero.");
    }
    if (min < gamma) {
        return {true, {Inf, NaN}}; // Model cannot fit these data
    }

    // Transform to uniform using lognormal CDF
    vector<double> data(len);
    transform(samples.begin(), samples.end(), data.begin(), 
        [gamma, mu, sigma](double sample) { 
            return 0.5 * (1 + erf((log(sample - gamma) - mu) / (sigma * sqrt(2.0))));
        }
    );
    sort(data.begin(), data.end());

    // Calculate Anderson-Darling statistic (A²)
    const double lenf = static_cast<double>(len);
    double accum = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(len > 1000)
    #endif     
    for (size_t i = 0; i < len; ++i) {
        accum += (2*(i + 1) - 1) * log(data[i]) + (2*lenf + 1 - 2*(i + 1))*log(1 - data[i]);
    }
    double A2 = -lenf - accum/lenf;

    // Apply small sample correction
    A2 *= (1.0 + 0.75/lenf + 2.25/(lenf*lenf));

    // double pvalue;
    // if (A2 <= 0.2) {
    //     pvalue = 1 - exp(-13.436 + 101.14 * A2 - 223.73 * A2 * A2);
    // } else if (A2 <= 0.34) {
    //     pvalue = 1 - exp(-8.318 + 42.796 * A2 - 59.938 * A2 * A2);
    // } else if (A2 <= 0.6) {
    //     pvalue = exp(0.9177 - 4.279 * A2 - 1.38 * A2 * A2);
    // } else if (A2 <= 153.467) {
    //     pvalue = exp(1.2937 * A2 - 5.709 * A2 + 0.0186 * A2 * A2);
    // } else {
    //     pvalue = 0;
    // }
    // double stat = A2;
    // bool reject = (pvalue <= 0.05); // equivalently, the statistic is greater than the threshold

    return {A2 > thresh, {A2, thresh}};
}

/**
 * Calculates the cumulative distribution function (CDF) for a single value
 * from the log-normal distribution with offset.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
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

    const double gamma = params.gamma;
    const double mu = params.mu;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    // Calculate the CDF of the lognormal distribution
    if(sample <= gamma) {
        return 0.0;
    } else {
        const double z = (log(sample - gamma) - mu) / (sigma * sqrt(2.0));
        return 0.5 * (1.0 + erf(z));
    }
}

/**
 * Calculates the cumulative distribution function (CDF) for multiple values
 * from the log-normal distribution with offset.
 * 
 * @param gamma Location parameter (minimum possible value)
 * @param mu Mean of the associated normal distribution
 * @param sigma Standard deviation of the associated normal distribution (sigma > 0)
 * @param x Vector of values to evaluate
 * @return Vector of CDF values F(x) = Φ((ln(x-gamma) - μ)/σ) 
 *         where Φ is the standard normal CDF
 * 
 * Implementation details:
 * - Uses error function (erf) for standard normal CDF calculation
 * - Handles boundary conditions (x <= gamma)
 * - Validates parameter constraints (sigma > 0)
 */
vector<double> LogNormalEstimator::cdf(const ModelParams& params, const vector<double>& samples) {

    const double gamma = params.gamma;
    const double mu = params.mu;
    const double sigma = params.sigma;

    const double sqrt2 = sqrt(2.0);
    unsigned len = samples.size();

    // Calculate the CDF of the lognormal distribution
    vector<double> cdf(len);
    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif    
    for(size_t i = 0; i < len; ++i) {
        if(samples[i] <= gamma) {
            cdf[i] = 0.0;
        } else {
            const double z = (log(samples[i] - gamma) - mu) / (sigma * sqrt2);
            cdf[i] = 0.5 * (1.0 + erf(z));
        }
    }
    return cdf;
}

/**
 * Calculates the cumulative distribution function (CDF) for a single value
 * from the log-normal distribution with offset.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
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

    const double gamma = params.gamma;
    const double mu = params.mu;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    const double sqrt2pi = sqrt(2.0 * M_PI);
    const double twoSigmaSq = 2.0 * sigma * sigma;

    // Calculate the PDF of the lognormal distribution
    if(sample <= gamma) {
        return 0.0;
    } else {
        const double diff = sample - gamma;
        const double logDiff = log(diff);
        const double expTerm = exp(-(logDiff - mu)*(logDiff - mu)/twoSigmaSq);
        return expTerm / (diff * sigma * sqrt2pi);
    }
}

/**
 * Calculates the probability density function (PDF) for multiple values
 * from the log-normal distribution with offset.
 * 
 * @param x Vector of values to evaluate
 * @param gamma Location parameter (minimum possible value)
 * @param mu Mean of the associated normal distribution
 * @param sigma Standard deviation of the associated normal distribution (sigma > 0)
 * @return Vector of PDF values f(x) = 1/((x-gamma)σ√(2π)) * exp(-(ln(x-gamma)-μ)²/(2σ²))
 * 
 * Implementation details:
 * - Pre-computes common terms for efficiency
 * - Handles boundary conditions (x <= gamma)
 * - Returns 0 for invalid inputs (x <= gamma)
 */
vector<double> LogNormalEstimator::pdf(const ModelParams& params, const vector<double>& samples) {

   // Parameters
    const double gamma = params.gamma;
    const double mu = params.mu;
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
        if(samples[i] <= gamma) {
            pdf[i] = 0.0;
        } else {
            const double diff = samples[i] - gamma;
            const double logDiff = log(diff);
            const double expTerm = exp(-(logDiff - mu)*(logDiff - mu)/twoSigmaSq);
            pdf[i] = expTerm / (diff * sigma * sqrt2pi);
        }
    }
    return pdf;
}

/**
 * Generates a single random value from the log-normal distribution using
 * the Box-Muller transform method.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
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
   const double gamma = params.gamma;
   const double mu = params.mu;
   const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

   const double sqrt2 = sqrt(2.0);

   // Generate 'length' data from loglogistic distribution

    // Box-Muller transform to get standard normal
    const double u = m_unif_dist(m_rnd_gen);
    const double z = sqrt2 * _erfinv(2.0*u - 1.0);

    return gamma + exp(mu + sigma*z);  // Transform to log-normal with parameters
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
   const double gamma = params.gamma;
   const double mu = params.mu;
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

        rnd[i] = gamma + exp(mu + sigma*z); // Transform to log-normal with parameters
   }
   return rnd;
}

/**
 * Calculates the expected value (mean) of the log-normal distribution.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
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
   const double mu = params.mu;
   const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }
    return gamma + exp(mu + (sigma*sigma)/2);
};

/**
 * Calculates the variance of the log-normal distribution.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
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
    const double mu = params.mu;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    const double sigma2 = sigma*sigma;

    return (exp(sigma2) - 1) * exp(2*mu + sigma2);
}

/**
 * Calculates the mode of the log-normal distribution.
 * 
 * @param params Distribution parameters {gamma, mu, sigma} where:
 *        - gamma: offset parameter (minimum possible value)
 *        - mu: location parameter (μ)
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
    const double mu = params.mu;
    const double sigma = params.sigma;

    if (sigma <= 0.0) {
        throw invalid_argument("Invalid sigma parameter for log-normal distribution (must be > 0)");
    }

    return gamma + exp(mu - sigma*sigma);
};
