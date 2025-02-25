#include "LogLogisticEstimator.h"
#include "alglib/optimization.h"

using namespace std;
using namespace alglib;

// LogLogisticEstimator class implementation
#define MAX_ITERS 100

struct optim_data {
    vector<double> samples;
    size_t len;
    double min;
};

/**
 * Calculates the function vector for log-logistic optimization
 * 
 * @param x Input parameter array [a, b, c] where:
 *          a: location parameter
 *          b: scale parameter
 *          c: shape parameter
 * @param fi Output function vector (3 components)
 * @param ptr Pointer to optimization parameters (samples and min value)
 */
void loglogistic_fvec(const real_1d_array& x, real_1d_array& fi, void* ptr) {

    // Get optimization parameters
    optim_data* p = (optim_data*)ptr;
    const vector<double>& samples = p->samples;
    const double min = p->min;
    const size_t len = p->len;

    // Get and validate parameters
    double a = max(eps, x[0]);
    double b = max(eps, x[1]);
    double c = max(eps, x[2]);
    a = std::min(a, min - eps);

    // Precompute common terms
    const double invc = 1.0/c;
    const double b2invc = pow(b, invc);
    const double logb = log(b);
    const double bc = b*c;
    const double nn = static_cast<double>(len);

    // Initialize sums
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

    // Main computation loop
    #pragma omp parallel for reduction(+:sum1,sum2,sum3) if(len > 1000)
    for (int i = 0; i < len; ++i) {
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
 * Calculates the Jacobian matrix for log-logistic optimization
 * 
 * @param x Input parameter array [a, b, c]
 * @param fi Output function vector (3 components)
 * @param jac Output Jacobian matrix (3x3)
 * @param ptr Pointer to optimization parameters
 */
void loglogistic_jac(const real_1d_array& x, real_1d_array& fi, real_2d_array& jac, void* ptr) {

    // Get optimization parameters
    optim_data* p = (optim_data*)ptr;
    const vector<double>& samples = p->samples;
    const double min = p->min;
    const size_t len = p->len;

    // Validate and initialize parameters
    double a = max(eps, x[0]);
    double b = max(eps, x[1]);
    double c = max(eps, x[2]);
    a = std::min(a, min - eps);

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
    #pragma omp parallel for reduction(+:sum1,sum2,sum3,f1a,f1b,f1c,f2a,f2b,f2c,f3a,f3b,f3c) if(len > 1000)
    for (int i = 0; i < len; ++i) {
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
        sum3 += logxma - 2.0*logxmab/aux1;

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

    // Final calculations for Jacobian elements
    f2a = -2.0*b2invc/(bc*c)*f2a;
    f2b = -nn/(bc*b) + 2.0*b2invc/(bc*bc)*f2b;
    f2c = -nn/(bc*c) + 2.0*b2invc/(bc*c*c)*f2c;
    f3a = invc*invc*f3a;
    f3b = -nn/(bc*c) - 2.0/(bc*c*c)*f3b;
    f3c = 2.0*nn*logb/(c*c*c) + nn/(c*c) - 2.0/(c*c*c)*f3c;

    // Calculate final function values
    sum2 = (nn - 2.0*b2invc*sum2)/(bc);
    sum3 = (-nn*(logb + c) + sum3)/(c*c);

    // Store function values
    fi[0] = sum1;
    fi[1] = sum2;
    fi[2] = sum3;

    // Store Jacobian matrix
    jac[0][0] = f1a; jac[0][1] = f1b; jac[0][2] = f1c;
    jac[1][0] = f2a; jac[1][1] = f2b; jac[1][2] = f2c;
    jac[2][0] = f3a; jac[2][1] = f3b; jac[2][2] = f3c;
}

LogLogisticEstimator::LogLogisticEstimator() : m_min_len(10) {
// Constructor
}

Model LogLogisticEstimator::fit(const vector<double>& samples) {

    size_t len = samples.size();

    vector<double> data = samples;
    sort(data.begin(), data.end());

    auto min = data[0];
    auto max = data[len-1];
    auto q1 = data[len/4];
    auto q3 = data[3*len/4];

    // sanity check
    if (len < m_min_len) {
        throw invalid_argument("Cannot fit anything with less than 'm_min_len' values.");
    }

    if (min <= 0) {
        throw invalid_argument("Cannot fit anything with min <= 0.");
    }

    // initial guess for a
    auto a0 = min - (max - min) / 1e4; // just a crude estimate of where the a could be, to start with
    if (a0 < 0) a0 = 0;

    // initial guess for b
    auto b0 = (max - min) / 2.0;
    if (b0 < 1e-6) b0 = 1e-6;  // to avoid some numerical errors in the fitting pass

    // initial guess for c using the 3rd ad 1st quartiles' relationship
    auto c0 = (q3 - q1) / (max - min);

    real_1d_array x;
    double x1[] = {a0, b0, c0};
    x.attach_to_ptr(3, x1);

    // Lower bounds
    real_1d_array bndl;
    double bndl1[] = {eps, eps, 0.05};
    bndl.attach_to_ptr(3, bndl1);

    // Upper bounds
    real_1d_array bndu;
    double bndu1[] = {min - eps, Inf, 0.5 - eps};
    bndu.attach_to_ptr(3, bndu1);

    real_1d_array s;
    double scales[] = {1e-6, 1e-6, 1e-6};
    s.attach_to_ptr(3, scales);

    // Optmization setup
    minlmstate state;
    minlmcreatevj(3, x, state);
    minlmsetbc(state, bndl, bndu);
    minlmsetcond(state, eps, MAX_ITERS);
    minlmsetscale(state, s);

    // Setup of OptGuard in debug only
    #ifdef _DEBUG
        minlmoptguardgradient(state, 0.001);
    #endif

    optim_data odata = {samples, len, min};

    // Optimization
    minlmoptimize(state, loglogistic_fvec, loglogistic_jac, NULL, &odata);

    // Results
    minlmreport rep;
    minlmresults(state, x, rep);

    if(rep.terminationtype < 0) {
        cerr << "Error: Optimization did not converge. Code: " << rep.terminationtype << endl;
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN}, {Inf, NAN}}
    }

    ModelParams params;
    params.a = x[0];
    params.b = x[1];
    params.c = x[2];

    auto [reject, gof] = this->gof(params, samples);
    if (!reject) {
        return {true, ModelType::LL3, params, gof};
    }
    else { 
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
    }
}

tuple<bool, GoF> LogLogisticEstimator::gof(const ModelParams& params, const vector<double>& samples) { // bool previous_model
// Anderson-Darling test the goodness of fit of the 3-loglogistic (A,B,C) for 
// the sample XS, from which the very parameters have been deduced, using a 0.05
// significance level, according to "Goodnes-of-fit techniques", D'Agostino and 
// Stephens (1986)

    bool previous_model = false;
    double a = params.a;
    double b = params.b;
    double c = params.c;
    double mu = std::log(b); // alpha in the book
    double thresh;
    if (previous_model) {
        thresh = 2.492; // for the case that parameters do not come from sample; 0.05, p. 105, table 4.2, right tail
    } else {
        thresh = 0.660; // threshold for the case of parameters coming from sample; 0.05 significance level; page 157, table 4.22, case 3
    }
    int len = samples.size();

    // sanity check
    double min = *min_element(samples.begin(), samples.end());
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
    if (!previous_model) {
        A2 *= (1.0 + 0.25 / len); // correction needed because both parameters are deduced from the same sample, table 4.22 case 3
    }

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
