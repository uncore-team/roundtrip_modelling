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
    //#pragma message("Compiling LogNormalEstimator with OpenMP support.")
#else
    #pragma message("Compiling LogNormalEstimator without OpenMP support.")
#endif
#include "LogNormalEstimator.h"
#include "alglib/optimization.h"

using namespace std;
using namespace alglib;

/**
 * Data structure for gamma parameter optimization using Modified Maximum Likelihood.
 * Contains references to input data and optimization parameters used in the
 * single-parameter optimization process.
 */
struct gamma_optim_data {
    const vector<double>& samples;
    const size_t& len;
    const double& min;
    const double& kr;
};

/**
 * Constructor for LogNormalEstimator class.
 * Initializes the base Estimator class with a minimum required sample size of 10.
 * This minimum size ensures reliable parameter estimation for the log-logistic distribution.
 */
LogNormalEstimator::LogNormalEstimator() : Estimator(10) {
}

/**
 * Evaluates objective function for gamma parameter optimization.
 * Implements Welford's online algorithm for stable mean and standard deviation calculation.
 * Assumes r=1 for the Modified Maximum Likelihood Estimation (MMLE-I).
 * 
 * @param samples Input vector of observations
 * @param min Minimum value in samples
 * @param kr Normal quantile for r/(n+1)
 * @param x Current gamma value being evaluated
 * @return Function value at x, or NaN if x is invalid
 * 
 * Implementation details:
 * - Uses Welford's algorithm for numerical stability
 * - Parallelizes computation for large datasets (>1000 samples)
 */
double LogNormalEstimator::gamma_fvec(const double gamma, const void* ptr) {

    // Get optimization parameters from the data structure
    gamma_optim_data* p = (gamma_optim_data*)ptr;
    const vector<double>& samples = p->samples;
    const size_t& len = p->len;
    const double& min = p->min;
    const double& kr = p->kr;

    // Validate input
    if (gamma >= min) {
        return NaN;
    }

    // Welford's online algorithm for mean and standard deviation
    double mu = 0.0;
    double M2 = 0.0;     // Second moment (used for standard deviation)
    double n = 0;        // Current count

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:mu,M2,n) if(len > OMP_THRESH)
    #endif
    for (const double& sample : samples) {
        const double value = log(sample - gamma);
        n++;
        const double delta = value - mu;  // Distance to current mean
        mu += delta/n;                    // Update mean incrementally
        M2 += delta*(value - mu);         // Update M2 incrementally
    }

    const double sigma = sqrt(M2/n); // Calculate standard deviation

    return log(min - gamma) - mu - kr*sigma;
}

/**
 * Fits a three-parameter log-normal distribution to sample data.
 * Uses a two-stage approach: first estimates gamma, then mu and sigma.
 * 
 * @param samples Vector of observations to fit
 * @return Model structure containing:
 *         - Validity flag
 *         - Distribution type (LN3)
 *         - Parameters {gamma, mu, sigma}
 *         - Goodness of fit metrics
 * 
 * Implementation details:
 * - First estimates offset parameter gamma using MMLE-I
 * - Then calculates mu and sigma using Welford's algorithm
 * - Validates fit using Anderson-Darling test
 * - Returns empty model if optimization fails or fit is poor
 */
Model LogNormalEstimator::fit(const vector<double>& samples) {

    const size_t len = samples.size();
    const double min = _min(samples);
    const double kr = sqrt(2.0)*_erfinv(2.0/(len + 1.0) - 1.0); // we assume r=1

    // initial guess for gamma
    double gamma = 0.5*min;

    // Setup search process
    int termcode; // Optimization exit code
    const double tolx = eps; // 
    const double a = 0;
    const double b = min - 1e-9;

    gamma_optim_data data = {samples, len, min, kr};
    tie(gamma, termcode) = fzero(gamma_fvec, a, b, tolx, &data);

    if(termcode < 0) {
        cerr << "Error: Optimization did not converge. Code: " << termcode << endl;
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
    }        

    // Calculate mean and standard deviation using Welford's online algorithm
    double mu = 0.0;
    double M2 = 0.0;
    double n = 0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:mu,M2,n) if(len > OMP_THRESH)
    #endif
    for (const double& sample : samples) {
        const double value = log(sample - gamma);
        n++;
        const double delta = value - mu;
        mu += delta/n;
        M2 += delta*(value - mu);
    }

    ModelParams params;
    params.gamma = gamma;
    params.mu = mu;
    params.sigma = sqrt(M2/(n - 1));

    auto [reject, gof_] = gof(params, samples);

    // Return model only if fit is acceptable
    // empty model = {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
    return reject ? Model() : Model{true, ModelType::LN3, params, gof_};
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
    const double sqrt2 = sqrt(2.0);
    transform(samples.begin(), samples.end(), data.begin(), 
        [gamma, mu, sigma, sqrt2](double sample) {
            const double w = (log(sample - gamma) - mu)/sigma;
            return 0.5*erfc(-w/sqrt2);
        }
    );
    sort(data.begin(), data.end());

    // Calculate Anderson-Darling statistic (A²)
    const double lenf = static_cast<double>(len);
    double accum = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(len > OMP_THRESH)
    #endif     
    for (size_t i = 0; i < len; ++i) {
        accum += (2*(i + 1) - 1) * log(data[i]) + (2*lenf + 1 - 2*(i + 1))*log(1 - data[i]);
    }
    double A2 = -lenf - accum/lenf;

    // Apply small sample correction
    A2 *= (1.0 + 0.75/lenf + 2.25/(lenf*lenf));

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
        #pragma omp parallel for if(len > OMP_THRESH)
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
        #pragma omp parallel for if(len > OMP_THRESH)
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
       #pragma omp parallel for if(length > OMP_THRESH)
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
}

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
}

/**
 * Finds a zero of a continuous function using a combination of bisection and interpolation.
 * Implementation based on the algorithm used in MATLAB's fzero function.
 * 
 * @param fun Function object to find zero for
 * @param a Left endpoint of initial interval
 * @param b Right endpoint of initial interval
 * @param tol Relative tolerance for convergence (default: 1e-6)
 * @param ptr Optional pointer to additional data needed by fun
 * 
 * @return pair<double, int> containing:
 *         - double: Approximate zero of the function
 *         - int: Exit flag indicating success/failure:
 *           1: Zero found
 *           -1: Terminated by output function
 *           -3: NaN or Inf encountered
 *           -4: Complex value encountered
 *           -5: Converged to singular point
 *           -6: No sign change detected in interval
 * 
 * Implementation details:
 * - Uses combination of bisection and inverse quadratic interpolation
 * - Ensures robust convergence through interval updating
 * - Handles corner cases (exact zeros, NaN/Inf values)
 * - Implements safeguards against numerical instability
 * 
 * @throws invalid_argument if a >= b (invalid interval)
 */
tuple<double, int> LogNormalEstimator::fzero(const function<double(const double, const void*)>& fun, double a, double b, double tol, void* ptr) {

    // Exit flags for different scenarios
    enum ExitFlag {
        FoundZero = 1,
        TerminatedByOutput = -1,
        NaNOrInfEncountered = -3,
        ComplexValueEncountered = -4,
        ConvergedToSingularPoint = -5,
        NoSignChangeDetected = -6
    };

    // Validate interval bounds
    if (a >= b) {
        throw invalid_argument("Left endpoint 'a' must be less than right endpoint 'b'");
    }

    // Evaluate function at interval endpoints
    double fa = fun(a, ptr);
    double fb = fun(b, ptr);

    // Check for exact solutions or invalid interval
    if (fa * fb >= 0) {
        if (fa == 0.0) return {a, ExitFlag::FoundZero};
        if (fb == 0.0) return {b, ExitFlag::FoundZero};
        return {0.0, ExitFlag::NoSignChangeDetected};
    }

    // Initialize algorithm variables
    double c = a;           // Previous step point
    double fc = fa;         // f(c)
    double d = b - a;       // Step size
    double e = d;           // Previous step size
    const int max_iter = 100;
    int iter = 0;

    // Main iteration loop
    while (fb != 0.0 && a != b) {
        // Update interval to ensure b is best point
        if (signbit(fb) == signbit(fc)) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }

        // Ensure b has smallest absolute value
        if (abs(fc) < abs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        // Check convergence
        const double m = 0.5 * (c - b);
        const double toler = 2.0 * tol * max(abs(b), 1.0);
        if (abs(m) <= toler || fb == 0.0) {
            break;
        }

        // Choose between bisection and interpolation
        if (abs(e) < toler || abs(fa) <= abs(fb)) {
            // Use bisection
            d = e = m;
        } else {
            // Use interpolation
            double s = fb / fa;
            double p, q;
            
            if (a == c) {
                // Linear interpolation
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                // Inverse quadratic interpolation
                const double r = fb / fc;
                const double s = fb / fa;
                p = s * (2.0 * m * (s - r) - (b - a) * (r - 1.0));
                q = (s - 1.0) * (r - 1.0) * (s - 1.0);
            }

            // Adjust sign
            if (p > 0.0) q = -q;
            else p = -p;

            // Accept interpolation if in bounds
            if (2.0 * p < 3.0 * m * q - abs(toler * q) && 
                p < abs(0.5 * e * q)) {
                e = d;
                d = p / q;
            } else {
                d = e = m;  // Use bisection
            }
        }

        // Update points for next iteration
        a = b;
        fa = fb;
        
        // Calculate new point with bounds checking
        if (abs(d) > toler) {
            b += d;
        } else {
            b += copysign(toler, m);
        }

        // Evaluate function at new point
        fb = fun(b, ptr);

        // Check for invalid results
        if (!isfinite(fb)) {
            return {b, isnan(fb) ? 
                ExitFlag::NaNOrInfEncountered : 
                ExitFlag::ComplexValueEncountered};
        }

        // Check iteration limit
        if (++iter >= max_iter) {
            return {b, ExitFlag::ConvergedToSingularPoint};
        }
    }

    // Return result
    return {b, ExitFlag::FoundZero};
}