/**
 * @brief Implementation of the LogLogisticEstimator class for three-parameter log-logistic distribution.
 * 
 * Provides methods for:
 * - Two-stage parameter estimation using maximum likelihood
 * - Shape parameter optimization with fixed location and scale
 * - Full three-parameter optimization using Levenberg-Marquardt
 * - Goodness of fit testing using Anderson-Darling
 * - OpenMP parallelization for large datasets
 * 
 * Based on D'Agostino & Stephens (1986), Chapter 4.
 */
#ifdef _OPENMP
    #include <omp.h>
    #pragma message("Compiling LogLogisticEstimator with OpenMP support.")
#else
    #pragma message("Compiling LogLogisticEstimator without OpenMP support.")
#endif
#include "LogLogisticEstimator.h"
#include "alglib/optimization.h"

using namespace std;
using namespace alglib;

/**
 * @brief Maximum number of iterations for optimization algorithms
 * 
 * Used in:
 * - c_estimate(): MAX_ITERS/2 iterations
 * - abc_estimate(): MAX_ITERS iterations
 * 
 * A larger value may improve convergence but increase computation time.
 */
#define MAX_ITERS 100

/**
 * Data structure for single-parameter optimization of shape parameter (c)
 * Used in c_estimate() and c_fvec() functions
 */
struct c_optim_data {
    double a;      // Location parameter (fixed during optimization)
    double b;      // Scale parameter (fixed during optimization)
    double mean;   // Target mean value from sample data
};

/**
 * Data structure for three-parameter optimization (a,b,c)
 * Used in abc_estimate(), abc_fvec() and abc_jac() functions
 */
struct abc_optim_data {
    vector<double> samples;  // Input data samples
    unsigned len;           // Number of samples
    double min;            // Minimum value in samples (used for parameter bounds)
};

/**
 * Optimized function for finding the shape parameter (c) of a log-logistic distribution
 * using a single-parameter optimization approach.
 * 
 * @param x Input parameter array containing a single value:
 *          x[0]: shape parameter (c) to be optimized
 * @param fi Output function vector (1 component) containing the objective function value
 * @param ptr Pointer to optimization data containing fixed parameters (a, b) and target mean
 * 
 * The function implements the following formula:
 * fi[0] = a + b*exp(lgamma(1+x) + lgamma(1-x) - lgamma(2)) - mean
 * 
 * where:
 * - a is the location parameter (fixed)
 * - b is the scale parameter (fixed)
 * - x is the shape parameter (to be optimized)
 * - mean is the target mean value
 * - lgamma is the log gamma function
 */
void c_fvec(const real_1d_array& x, real_1d_array& fi, void *ptr) {
    // Get optimization parameters from the data structure
    c_optim_data* p = (c_optim_data*)ptr;
    const double a = p->a;      // Location parameter (fixed)
    const double b = p->b;      // Scale parameter (fixed)
    const double mean = p->mean;// Target mean value

    // Calculate terms for the gamma functions
    const double x1 = 1.0 + x[0];  // For gamma(1+x)
    const double x2 = 1.0 - x[0];  // For gamma(1-x)

    // Calculate objective function:
    // The difference between the theoretical mean and the sample mean
    fi[0] = a + b*exp(lgamma(x1) + lgamma(x2) - lgamma(x1 + x2)) - mean;
}

/**
 * Calculates the function vector for log-logistic optimization using maximum likelihood estimation.
 * This function computes the gradient of the log-likelihood function with respect to the parameters.
 * 
 * @param x Input parameter array [a, b, c] where:
 *          a: location parameter (minimum possible value)
 *          b: scale parameter (controls spread)
 *          c: shape parameter (controls skewness)
 * @param fi Output function vector (3 components) containing partial derivatives:
 *          fi[0]: derivative with respect to a
 *          fi[1]: derivative with respect to b
 *          fi[2]: derivative with respect to c
 * @param ptr Pointer to optimization data structure containing:
 *          - samples: vector of observed values
 *          - len: number of samples
 *          - min: minimum value in samples
 * 
 * Implementation details:
 * - Uses OpenMP for parallel computation when dataset is large (>1000 samples)
 * - Parameters are bounded to ensure numerical stability
 * - Precomputes common terms to improve performance
 * - Uses reduction for parallel accumulation of sums
 * 
 * The function implements the maximum likelihood equations:
 * ∂L/∂a = sum(terms involving a)
 * ∂L/∂b = sum(terms involving b)
 * ∂L/∂c = sum(terms involving c)
 * 
 * where L is the log-likelihood function for the log-logistic distribution.
 */
void abc_fvec(const real_1d_array& x, real_1d_array& fi, void* ptr) {

    // Get optimization parameters
    abc_optim_data* p = (abc_optim_data*)ptr;
    const vector<double>& samples = p->samples;
    const double min = p->min;
    const unsigned len = p->len;

    // Get and validate parameters
    double a = MIN( MAX(eps, x[0]), min - eps);
    double b = MAX(eps, x[1]);
    double c = MAX(eps, x[2]);

    // Precompute common terms
    const double invc = 1.0/c;
    const double b2invc = pow(b, invc);
    const double logb = log(b);
    const double bc = b*c;
    const double nn = static_cast<double>(len);

    // Initialize sums
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

    // Main computation loop
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:sum1,sum2,sum3) if(len > 1000)
    #endif
    for (unsigned i = 0; i < len; ++i) {
        const double xma = samples[i] - a;	
        const double xma2invc = pow(xma, invc);
        const double logxma = log(xma);
        const double xmab = xma/b;
        const double aux = xma2invc + b2invc;
        const double logxmab = log(xmab);

        // Update partial sums
		sum1 += ((1.0 + invc) - (2.0*b2invc/c)/aux)/xma;
		sum2 += 1.0/aux;
		sum3 += logxma - 2.0*logxmab/(pow(xmab, b2invc) + 1.0);
	}

    // Final calculations
	sum2 = (nn - 2.0*b2invc*sum2)/(bc);
	sum3 = (-nn*(logb + c) + sum3)/(c*c);

    // Store results
    fi[0] = sum1;
    fi[1] = sum2;
    fi[2] = sum3;
}

/**
 * Calculates the Jacobian matrix for log-logistic optimization using maximum likelihood estimation.
 * This function computes the second-order derivatives (Hessian matrix) of the log-likelihood function.
 * 
 * @param x Input parameter array [a, b, c] where:
 *          a: location parameter (minimum possible value)
 *          b: scale parameter (controls spread)
 *          c: shape parameter (controls skewness)
 * @param fi Output function vector (3 components) containing first derivatives
 * @param jac Output Jacobian matrix (3x3) containing second derivatives:
 *          [∂²L/∂a², ∂²L/∂a∂b, ∂²L/∂a∂c]
 *          [∂²L/∂b∂a, ∂²L/∂b², ∂²L/∂b∂c]
 *          [∂²L/∂c∂a, ∂²L/∂c∂b, ∂²L/∂c²]
 * @param ptr Pointer to optimization data structure containing:
 *          - samples: vector of observed values
 *          - len: number of samples
 *          - min: minimum value in samples
 * 
 * Implementation details:
 * - Uses OpenMP for parallel computation when dataset is large (>1000 samples)
 * - Parameters are bounded to ensure numerical stability
 * - Precomputes common terms to improve performance
 * - Uses reduction for parallel accumulation of sums and matrix elements
 * - Implements analytical derivatives for better numerical accuracy
 * 
 * The function calculates both:
 * 1. First derivatives (stored in fi vector)
 * 2. Second derivatives (stored in jac matrix)
 * 
 * These are used by the Levenberg-Marquardt optimizer to determine:
 * - Direction of steepest descent
 * - Step size for parameter updates
 * - Convergence criteria
 */
void abc_jac(const real_1d_array& x, real_1d_array& fi, real_2d_array& jac, void* ptr) {

    // Get optimization parameters
    abc_optim_data* p = (abc_optim_data*)ptr;
    const vector<double>& samples = p->samples;
    const double min = p->min;
    const unsigned len = p->len;

    // Validate and initialize parameters
    double a = MIN(MAX(eps, x[0]), min - eps);
    double b = MAX(eps, x[1]);
    double c = MAX(eps, x[2]);

    // Precompute common terms
    const double invc = 1.0/c;
    const double b2invc = pow(b, invc);
    const double logb = log(b);
    const double bc = b*c;
    const double nn = static_cast<double>(len);
    const double c2 = c*c;
    const double c3 = c2*c;

    // Initialize accumulators for sums and Jacobian elements
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    double f1a = 0.0, f1b = 0.0, f1c = 0.0;
    double f2a = 0.0, f2b = 0.0, f2c = 0.0;
    double f3a = 0.0, f3b = 0.0, f3c = 0.0;

    // Main computation loop with OpenMP parallelization for large datasets
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:sum1,sum2,sum3,f1a,f1b,f1c,f2a,f2b,f2c,f3a,f3b,f3c) if(len > 1000)
    #endif
    for (unsigned i = 0; i < len; ++i) {
        // Compute intermediate values
        const double xma = samples[i] - a;
        const double xma2invc = pow(xma, invc);
        const double logxma = log(xma);
        const double xma2 = xma * xma;
        const double xmab = xma/b;
        const double xmab2invc = pow(xmab, invc);
        const double logxmab = log(xmab);

        // Compute auxiliary terms
        const double aux = xma2invc + b2invc;
        const double aux2 = aux * aux;
        const double aux1 = xmab2invc + 1.0;
        const double aux12 = aux1 * aux1;

        // Update sums for function values
        sum1 += ((1.0 + invc) - (2.0*b2invc/c)/aux)/xma;
        sum2 += 1.0/aux;
        sum3 += logxma - 2.0*logxmab/(pow(xmab, b2invc) + 1.0);

        // Update Jacobian elements
        f1a += (1.0 + invc)/xma2 - 2.0*b2invc/c*(c*aux + xma2invc)/(c*xma2*aux2);
        f1b += -2.0/(c*xma)*b2invc*xma2invc/(bc*aux2);
        f1c += -1.0/(c2*xma) + 2.0*b2invc*(c*aux + logb*xma2invc - xma2invc*logxma)/(c3*xma*aux2);
        
        f2a += xma2invc/(xma*aux2);
        f2b += (c*aux - xma2invc)/aux2;
        f2c += (c*aux + logb*xma2invc - xma2invc*logxma)/aux2;
        
        f3a += -1.0/xma + 2.0*(c*xmab2invc - xmab2invc*logxmab + c)/(c*xma*aux12);
        f3b += (xmab2invc*logxmab - c*aux1)/aux12;
        f3c += logxma + (logxmab*(xmab2invc*logxmab - 2.0*c*aux1))/(c*aux12);
    }

    // Calculate final function values
    sum2 = (nn - 2.0*b2invc*sum2)/(bc);
    sum3 = (-nn*(logb + c) + sum3)/c2;

    // Final calculations for Jacobian elements
    f2a = -2.0*b2invc/(bc*c)*f2a;
    f2b = -nn/(bc*b) + 2.0*b2invc/(bc*bc)*f2b;
    f2c = -nn/(bc*c) + 2.0*b2invc/(bc*c2)*f2c;
    f3a = invc*invc*f3a;
    f3b = -nn/(bc*c) - 2.0/(bc*c2)*f3b;
    f3c = 2.0*nn*logb/c3 + nn/c2 - 2.0/c3*f3c;

    // Store function values
    fi[0] = sum1;
    fi[1] = sum2;
    fi[2] = sum3;

    // Store Jacobian matrix
    jac[0][0] = f1a; jac[0][1] = f1b; jac[0][2] = f1c;
    jac[1][0] = f2a; jac[1][1] = f2b; jac[1][2] = f2c;
    jac[2][0] = f3a; jac[2][1] = f3b; jac[2][2] = f3c;
}

/**
 * Constructor for LogLogisticEstimator class.
 * Initializes the base Estimator class with a minimum required sample size of 10.
 * This minimum size ensures reliable parameter estimation for the log-logistic distribution.
 */
LogLogisticEstimator::LogLogisticEstimator() : Estimator(10) {
}

/**
 * Estimates the shape parameter (c) of a log-logistic distribution while keeping location (a) and scale (b) fixed.
 * Uses Levenberg-Marquardt optimization with single-parameter bounded optimization.
 * 
 * @param samples Input vector of observed values
 * @param a Fixed location parameter (minimum possible value)
 * @param b Fixed scale parameter (controls spread)
 * @param c Input/Output shape parameter (controls skewness):
 *          - Input: Initial guess for optimization
 *          - Output: Optimized value if successful
 * 
 * @return Termination type from optimizer:
 *         >0: successful completion
 *         =0: maximum number of iterations reached
 *         <0: optimization failed
 * 
 * Implementation details:
 * - Uses bounded optimization to keep c in [0.05, 0.5-eps]
 * - Scales parameter for numerical stability
 * - Uses mean-matching approach for optimization
 * - Maximum iterations set to MAX_ITERS/2
 */
int LogLogisticEstimator::c_estimate(const vector<double>& samples, const double& a, const double& b,  double& c) {

    //unsigned len = samples.size();
    double mean = _mean(samples);
    c_optim_data data = {a, b, mean};

    // Initial guess of c
    real_1d_array x;
    double x1[] = {c};
    x.attach_to_ptr(1, x1);

    // Lower/Upper bounds
    real_1d_array bndl, bndu;
    double bndl1[] = {0.05};
    bndl.attach_to_ptr(1, bndl1);
    double bndu1[] = {0.5 - eps};
    bndu.attach_to_ptr(1, bndu1);

    // Optmization setup
    minlmstate state;
    minlmcreatev(1, x, 1e-3, state);
    minlmsetbc(state, bndl, bndu);
    minlmsetcond(state, eps, MAX_ITERS / 2);

    real_1d_array s;
    double scales[] = {1e-6};
    s.attach_to_ptr(1, scales);
    minlmsetscale(state, s);

    // Optimization
    minlmoptimize(state, c_fvec, NULL, &data);

    // Results
    minlmreport rep;
    minlmresults(state, x, rep);
    c = x[0];

    return rep.terminationtype;
}

/**
 * Estimates all three parameters (a,b,c) of a log-logistic distribution simultaneously.
 * Uses Levenberg-Marquardt optimization with multi-parameter bounded optimization.
 * 
 * @param samples Input vector of observed values
 * @param a Input/Output location parameter:
 *          - Input: Initial guess
 *          - Output: Optimized value
 * @param b Input/Output scale parameter:
 *          - Input: Initial guess
 *          - Output: Optimized value
 * @param c Input/Output shape parameter:
 *          - Input: Initial guess
 *          - Output: Optimized value
 * 
 * @return Termination type from optimizer:
 *         >0: successful completion
 *         =0: maximum number of iterations reached
 *         <0: optimization failed
 * 
 * Implementation details:
 * - Uses bounded optimization:
 *   * a ∈ [eps, min-eps]
 *   * b ∈ [eps, Inf]
 *   * c ∈ [0.05, 0.5-eps]
 * - Uses both function values and Jacobian matrix
 * - Implements maximum likelihood estimation
 * - Maximum iterations set to MAX_ITERS
 * - Uses parameter scaling for numerical stability
 */
int LogLogisticEstimator::abc_estimate(const vector<double>& samples, double& a, double& b,  double& c) {

    unsigned len = samples.size();
    double min = _min(samples);
    abc_optim_data data = {samples, len, min};

    // Initial guess of {a, b, c}
    real_1d_array x;
    double x1[] = {a, b, c};
    x.attach_to_ptr(3, x1);

    // Lower/Upper bounds
    real_1d_array bndl, bndu;
    double bndl1[] = {eps, eps, 0.05};
    bndl.attach_to_ptr(3, bndl1);
    double bndu1[] = {min - eps, Inf, 0.5 - eps};
    bndu.attach_to_ptr(3, bndu1);

    // Optmization setup
    minlmstate state;
    minlmcreatevj(3, x, state);
    minlmsetbc(state, bndl, bndu);
    minlmsetcond(state, eps, MAX_ITERS);

    real_1d_array s;
    double scales[] = {1e-6, 1e-6, 1e-6};
    s.attach_to_ptr(3, scales);
    minlmsetscale(state, s);

    // Optmization
    minlmoptimize(state, abc_fvec, abc_jac, NULL, &data);

    // Results
    minlmreport rep;
    minlmresults(state, x, rep);
    a = x[0];
    b = x[1];
    c = x[2];

    return rep.terminationtype;
}

/**
 * Fits a three-parameter log-logistic distribution to the given samples.
 * Uses a two-stage optimization approach for better convergence.
 * 
 * @param samples Input vector of observed values (must be positive)
 * 
 * @return Model structure containing:
 *         - defined: true if fit was successful
 *         - type: ModelType::LL3 for successful fit, ModelType::None otherwise
 *         - params: {a, b, c} parameters of the distribution
 *         - gof: Goodness of fit statistics
 * 
 * Implementation details:
 * - First stage: Estimates shape parameter c while keeping a and b fixed
 * - Second stage: Refines all parameters (a, b, c) simultaneously
 * - Initial guesses:
 *   * a: slightly below minimum sample value
 *   * b: median of (samples - a)
 *   * c: 0.25 (midpoint of typical range)
 * - Performs goodness of fit test after parameter estimation
 * 
 * @throws invalid_argument if:
 *         - samples.size() < m_min_len
 *         - min(samples) <= 0
 */
Model LogLogisticEstimator::fit(const vector<double>& samples) {

    unsigned len = samples.size();
    double min = _min(samples);
    double max = _max(samples);

    // sanity check
    if (len < m_min_len) {
        throw invalid_argument("Cannot fit anything with less than 'm_min_len' values.");
    }
    if (min <= 0) {
        throw invalid_argument("Cannot fit anything with min <= 0.");
    }

    // initial guess for a
    double a = min - (max - min) / 1e4; // just a crude estimate of where the a could be, to start with
    a = _max(a, 0.0);

    // initial guess for b
    vector<double> data(len);
    transform(samples.begin(), samples.end(), data.begin(), [a](double value){ return value - a; });

    double b = _median(data);
    b = _max(b, 1e-6);  // to avoid some numerical errors in the fitting pass

    // initial guess for c (shape parameter)
    double c = 0.25; // midpoint between typical bounds (0.05, 0.5)

    // estimate c parameter first
    int termcode = c_estimate(samples, a, b, c);
    if(termcode < 0 ) {
        cerr << "Error: Optimization did not converge. Code: " << termcode << endl;
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
    }

    // estimate a, b, and c parameters together
    termcode = abc_estimate(samples, a, b, c);
    if(termcode < 0) {
        cerr << "Error: Optimization did not converge. Code: " << termcode << endl;
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
    }

    ModelParams params;
    params.a = a;
    params.b = b;
    params.c = c;

    auto [reject, gof_] = gof(params, samples);
    if (!reject) {
        return {true, ModelType::LL3, params, gof_};
    }
    else { 
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
    }
}

/**
 * Performs Anderson-Darling goodness of fit test for log-logistic distribution.
 * Based on "Goodness-of-fit techniques" by D'Agostino and Stephens (1986).
 * 
 * @param params Model parameters {a, b, c} to test
 * @param samples Input vector of observed values
 * 
 * @return tuple<bool, GoF> containing:
 *         - bool: true if null hypothesis should be rejected (poor fit)
 *         - GoF: {statistic, threshold} for the Anderson-Darling test
 * 
 * Implementation details:
 * - Uses Case 3 from Table 4.22 (p.157) for parameters estimated from same sample
 * - Significance level: 0.05 (95% confidence)
 * - Applies small sample size correction factor
 * - Transforms data to uniform distribution for testing
 * 
 * Rejection criteria:
 * - If any sample ≤ location parameter a
 * - If test statistic > threshold (0.660 for parameters from same sample)
 */
tuple<bool, GoF> LogLogisticEstimator::gof(const ModelParams& params, const vector<double>& samples) { // bool previous_model

    // bool previous_model = false;
    double a = params.a;
    double b = params.b;
    double c = params.c;
    double mu = log(b); // alpha in the book
    double thresh;
    // if (previous_model) {
    //     thresh = 2.492; // for the case that parameters do not come from sample; 0.05, p. 105, table 4.2, right tail
    // } else {
        thresh = 0.660; // threshold for the case of parameters coming from sample; 0.05 significance level; page 157, table 4.22, case 3
    // }
    int len = samples.size();

    // sanity check
    double min = _min(samples);
    if (min <= a) {
        return {true, {Inf, thresh}}; // cannot accept a distribution if some value falls below its minimum
    }

    // transform sample to LL2 (0 offset) and model (a,b,c) into C++ model (mu,sigma)
    vector<double> data(len);
    transform(
        samples.begin(), samples.end(),
        data.begin(),
        [a, c, mu](double sample) { return 1.0 / (1.0 + exp(-(log(sample - a) - mu) / c)); }
    );
    sort(data.begin(), data.end());

    // calculate statistic: A2 for case 3 (both parameters were deduced from
    // the same sample). This statistic measures the squared distance
    // between experimental and theoretical Zs, and, indirectly, between 
    // theoretical and experimental Xs (p. 100)
    double accum = 0.0;
    for (int i = 0; i < len; ++i) {
        accum += (2 * (i + 1) - 1) * log(data[i]) + (2 * len + 1 - 2 * (i + 1)) * log(1 - data[i]);
    }

    double A2 = -len - (1.0 / len) * accum;
    // if (!previous_model) {
        A2 *= (1.0 + 0.25 / len); // correction needed because both parameters are deduced from the same sample, table 4.22 case 3
    // }

    double stat = A2; // this statistic follows certain right-tailed distribution. We can set in
                      // that distribution a threshold value (in its support)
                      // corresponding to a given significance level. 
                      // Then, if the value calculated for the statistic falls 
                      // above that threshold, the hypothesis should be rejected
                      // (this is easier as the significance level grows).
                      // The p-value is the probability of the statistic distribution to
                      // produce a value of the statistic equal to or greater than the
                      // calculated one. The p-value will shrink as more strongly rejected
                      // is the null hypothesis. We do not calculate it here
                      // because the distribution of the statistic is not found in
                      // the book.

    // test the hypothesis 
    bool reject = (stat > thresh); // equivalently, the p-value is smaller than the significant level

    return {reject, {stat, thresh}};
}

/**
 * Calculates the cumulative distribution function (CDF) for a single value from 
 * the log-logistic distribution. Implements the analytical form of the CDF based
 * on the three-parameter log-logistic distribution.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a (offset): location parameter (minimum possible value)
 *        - b (alpha): scale parameter
 *        - c (sigma): shape parameter
 * @param sample Single input value to evaluate
 * 
 * @return CDF value F(x) = 1/(1 + ((x-a)/b)^(-1/c))
 * 
 * Implementation details:
 * - Handles boundary conditions (x <= a) explicitly
 * - Uses location-scale parameterization
 * - Returns 0 for values <= location parameter
 */
double LogLogisticEstimator::cdf(const ModelParams& params, const double& sample) {

    const double offset = params.a;
    const double alpha = params.b;
    const double beta = 1.0 / params.c;

    if (sample <= offset) {
        return 0.0;  // For x <= a, CDF is 0
    }
    return 1.0 / (1.0 + pow((sample - offset) / alpha, -beta));
}

/**
 * Calculates the cumulative distribution function (CDF) for multiple values from
 * the log-logistic distribution. Provides vectorized implementation for efficiency
 * on large datasets.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a (offset): location parameter (minimum possible value)
 *        - b (alpha): scale parameter
 *        - c (sigma): shape parameter
 * @param samples Vector of input values to evaluate
 * 
 * @return Vector of CDF values F(x) = 1/(1 + ((x-a)/b)^(-1/c))
 * 
 * Implementation details:
 * - Uses OpenMP for parallel computation when dataset is large (>1000 samples)
 * - Pre-computes common terms outside the loop for efficiency
 * - Handles boundary conditions (x <= a) explicitly
 * - Returns vector of same size as input samples
 */
vector<double> LogLogisticEstimator::cdf(const ModelParams& params, const vector<double>& samples) {

    const double offset = params.a;
    const double alpha = params.b;
    const double beta = 1.0 / params.c;

    size_t len = samples.size();
    vector<double> cdf(len);
    
    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif
    for (size_t i = 0; i < len; ++i) {
        if (samples[i] <= offset) {
            cdf[i] = 0.0;
        } else {
            cdf[i] = 1.0 / (1.0 + pow((samples[i] - offset) / alpha, -beta));
        }
    }
    return cdf;
}

/**
 * Calculates the probability density function (PDF) for a single value from
 * the log-logistic distribution. Implements the analytical form of the PDF.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a (offset): location parameter (minimum possible value)
 *        - b (alpha): scale parameter
 *        - c (sigma): shape parameter
 * @param sample Single input value to evaluate
 * 
 * @return PDF value f(x) = (β/α)((x-a)/α)^(β-1)/(1+((x-a)/α)^β)^2
 * where β = 1/c
 * 
 * Implementation details:
 * - Pre-computes common terms for efficiency
 * - Handles boundary conditions (x <= a) explicitly
 * - Uses transformed parameters for numerical stability
 */
double LogLogisticEstimator::pdf(const ModelParams& params, const double& sample) {
    const double offset = params.a;
    const double alpha = params.b;
    const double sigma = params.c;
    const double beta = 1.0 / params.c;

    if (sample <= offset) {
        return 0.0;
    }

    const double temp = pow((sample - offset) / alpha, -beta);
    return temp / (sigma * (sample - offset) * pow(1.0 + temp, 2));
}

/**
 * Calculates the probability density function (PDF) for multiple values from
 * the log-logistic distribution. Provides vectorized implementation for efficiency.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a (offset): location parameter (minimum possible value)
 *        - b (alpha): scale parameter
 *        - c (sigma): shape parameter
 * @param samples Vector of input values to evaluate
 * 
 * @return Vector of PDF values f(x) = (β/α)((x-a)/α)^(β-1)/(1+((x-a)/α)^β)^2
 * where β = 1/c
 * 
 * Implementation details:
 * - Uses OpenMP for parallel computation when dataset is large (>1000 samples)
 * - Pre-computes common terms outside the loop for efficiency
 * - Handles boundary conditions (x <= a) explicitly
 */
vector<double> LogLogisticEstimator::pdf(const ModelParams& params, const vector<double>& samples) {
    // Get and validate parameters
    const double offset = params.a;
    const double alpha = params.b;
    const double sigma = params.c;
    const double beta = 1.0 / params.c;

    // Calculate the PDF for the loglogistic distribution
    const size_t len = samples.size();
    vector<double> pdf(len);

    #ifdef _OPENMP
        #pragma omp parallel for if(len > 1000)
    #endif
    for (size_t i = 0; i < len; ++i) {
        if (samples[i] <= offset) {
            pdf[i] = 0;  // For x <= a, PDF is 0
        } else {
            const double temp = pow((samples[i] - offset) / alpha, -beta);
            pdf[i] = temp / (sigma * (samples[i] - offset) * pow(1.0 + temp, 2));
        }
    }
    return pdf;
}

/**
 * Generates a single random value from the log-logistic distribution using
 * the inverse transform sampling method. Uses the uniform distribution
 * inherited from the base class.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a (offset): location parameter (minimum possible value)
 *        - b (alpha): scale parameter
 *        - c (sigma): shape parameter
 * @return Random value following the log-logistic distribution
 * 
 * Uses inverse transform sampling with the formula:
 * X = a + exp(μ + σ*log(p/(1-p)))
 * where:
 * - μ = log(b)
 * - σ = c
 * - p ~ U(0,1)
 * 
 * Implementation details:
 * - Uses uniform distribution from base class
 * - Pre-computes common terms for efficiency
 * - Applies location-scale transformation
 */
double LogLogisticEstimator::rnd(const ModelParams& params) {
    const double offset = params.a;
    const double mu = log(params.b);
    const double sigma = params.c;

    const double p = m_unif_dist(m_rnd_gen);
    return offset + exp(mu + sigma * log(p / (1.0 - p)));
}

/**
 * Generates multiple random values from the log-logistic distribution using
 * vectorized inverse transform sampling. Optimized for generating large
 * numbers of random values efficiently.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a (offset): location parameter (minimum possible value)
 *        - b (alpha): scale parameter
 *        - c (sigma): shape parameter
 * @param length Number of random values to generate
 * @return Vector of random values following the log-logistic distribution
 * 
 * Uses inverse transform sampling with the formula:
 * X = a + exp(μ + σ*log(p/(1-p)))
 * where:
 * - μ = log(b)
 * - σ = c
 * - p ~ U(0,1)
 * 
 * Implementation details:
 * - Uses OpenMP for parallel generation when length > 1000
 * - Reuses uniform distribution from base class
 * - Pre-computes common terms for efficiency
 */
vector<double> LogLogisticEstimator::rnd(const ModelParams& params, const unsigned& length) {
    // Parameters
    const double offset = params.a;
    const double alpha = params.b;
    const double sigma = params.c;
    const double mu = log(alpha);

    // Generate 'length' data from loglogistic distribution
    vector<double> rnd(length);

    #ifdef _OPENMP
        #pragma omp parallel for if(length > 1000)
    #endif
    for (unsigned i = 0; i < length; ++i) {
        const double p = m_unif_dist(m_rnd_gen);  // Generate random number from uniform distribution
        rnd[i] = offset + exp(mu + sigma * log(p / (1.0 - p)));
    }
    return rnd;
}

/**
 * Calculates the expectation (mean) of the log-logistic distribution.
 * For shape parameter c > 1, the expectation exists and is given by
 * the analytical formula.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a: location parameter (minimum possible value)
 *        - b: scale parameter
 *        - c: shape parameter
 * @return Expected value E[X] = a + b*π*c/sin(π*c) when c < 1, NaN otherwise
 * 
 * Implementation details:
 * - Validates existence condition (c < 1)
 * - Uses transformed parameters for numerical stability
 * - Returns Inf when expectation does not exist
 */
double LogLogisticEstimator::expectation(const ModelParams& params) {
    const double c = params.c;
    const double pic = M_PI*c;
    if (c >= 1.0) {
        return NaN;  // Expectation does not exist for c >= 1
    }
    return params.a + params.b * pic / sin(pic);
}

/**
 * Calculates the variance of the log-logistic distribution.
 * For shape parameter c > 2, the variance exists and is given by
 * the analytical formula.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a: location parameter (minimum possible value)
 *        - b: scale parameter
 *        - c: shape parameter
 * @return Variance when c < 0.5, NaN otherwise
 * 
 * Implementation details:
 * - Validates existence condition (c < 0.5)
 * - Pre-computes common terms for efficiency
 * - Returns NaN when variance does not exist
 */
double LogLogisticEstimator::variance(const ModelParams& params) {

    const double c = params.c;
    if (c >= 0.5) {
        return NaN;  // Variance does not exist for c >= 0.5
    }

    const double pic = M_PI*c;
    return params.b * (2*pic/(sin(2*pic)) - pow(pic/sin(pic), 2.0));
}

/**
 * Calculates the mode of the log-logistic distribution.
 * The mode exists for all parameter values and has an analytical form.
 * 
 * @param params Distribution parameters {a, b, c} where:
 *        - a: location parameter (minimum possible value)
 *        - b: scale parameter
 *        - c: shape parameter
 * @return Mode = a + b*((1-c)/(1+c))^c when c < 1, a otherwise
 * 
 * Implementation details:
 * - Validates existence condition (c < 1)
 * - Uses transformed parameters for numerical stability
 */
double LogLogisticEstimator::mode(const ModelParams& params) {
    const double c = params.c;
    if (c >= 1) {
        return NaN;  // Mode does not exist for c >= 1
    }
    return params.a + params.b * pow((1-c)/(1+c), c);
}