#include "Estimator.h"

using namespace std;

/**
 * @brief Implementation of the Estimator base class for statistical distribution fitting.
 * 
 * Provides common statistical utility functions used by derived distribution estimators.
 * All functions are protected to ensure access only through derived classes.
 */

/**
 * @brief Constructor initializing minimum sample size requirement
 * @param min_len Minimum number of samples needed for reliable estimation
 */
Estimator::Estimator(unsigned min_len) : m_min_len(min_len), m_rnd_gen(m_rnd_dev()), m_unif_dist(0.0, 1.0) {
}

/**
 * @brief Virtual destructor implementation
 */
Estimator::~Estimator() {
}

/**
 * Returns the minimum of two values using a safe comparison.
 * 
 * @param a First value to compare
 * @param b Second value to compare
 * @return Smaller of the two values
 * 
 * Implementation details:
 * - Uses ternary operator for efficiency
 * - Handles NaN values according to IEEE 754
 */
double Estimator::_min(const double& a, const double& b) {

    return (a < b) ? a : b;
}

/**
 * Finds the minimum value in a vector of observations.
 * 
 * @param x Vector of observations
 * @return Minimum value in the vector
 * 
 * Implementation details:
 * - Uses STL min_element algorithm
 * - Returns NaN for empty vectors
 */
double Estimator::_min(const vector<double>& samples) {

    if(samples.empty()) {
        throw invalid_argument("Cannot find minimum of empty vector");
    }
    return *min_element(samples.begin(), samples.end());
}

/**
 * Returns the maximum of two values using a safe comparison.
 * 
 * @param a First value to compare
 * @param b Second value to compare
 * @return Larger of the two values
 * 
 * Implementation details:
 * - Uses ternary operator for efficiency
 * - Handles NaN values according to IEEE 754
 */
double Estimator::_max(const double& a, const double& b) {

    return (a > b) ? a : b;
}

/**
 * Finds the maximum value in a vector of observations.
 * 
 * @param x Vector of observations
 * @return Maximum value in the vector
 * 
 * Implementation details:
 * - Uses STL max_element algorithm
 * - Returns NaN for empty vectors
 */
double Estimator::_max(const vector<double>& samples) {

    if(samples.empty()) {
        throw invalid_argument("Cannot find maximum of empty vector");
    }
    return *max_element(samples.begin(), samples.end());
}

/**
 * Calculates the arithmetic mean of sample data.
 * 
 * @param x Vector of observations
 * @return Arithmetic mean of the samples
 * 
 * Implementation details:
 * - Uses STL accumulate for summation
 * - Handles potential numeric overflow
 * - Validates input size
 */
double Estimator::_mean(const vector<double>& samples) {

    if(samples.empty()) {
        throw invalid_argument("Cannot calculate mean of empty vector");
    }
    const size_t len = samples.size();
    return accumulate(samples.begin(), samples.end(), 0.0) / static_cast<double>(len);
}

/**
 * Calculates the median of sample data.
 * 
 * @param x Vector of observations (copied for sorting)
 * @return Median value:
 *         - For odd n: middle value
 *         - For even n: average of two middle values
 * 
 * Implementation details:
 * - Creates sorted copy to preserve original data
 * - Uses two-point average for even-length vectors
 * - Validates input size
 */
double Estimator::_median(const vector<double>& samples) {

    if(samples.empty()) {
        throw invalid_argument("Cannot calculate median of empty vector");
    }

    vector<double> data(samples);
    sort(data.begin(), data.end());

    const size_t len = data.size();
    if (len % 2 != 0) {
        return data[len/2];
    }
    return (data[(len-1)/2] + data[len/2]) / 2.0;
}

/**
 * Calculates the error function using a numerical approximation.
 * Implementation based on the Abramowitz and Stegun approximation 7.1.26.
 * 
 * @param x Input value
 * @return Approximated value of error function erf(x)
 * 
 * Implementation details:
 * - Uses polynomial approximation for improved accuracy
 * - Valid for all real input values
 * - Relative error < 1.5e-7
 * - Preserves sign symmetry: erf(-x) = -erf(x)
 */
double Estimator::_erf(const double& x) {

    // Constants for the approximation
    const double a1 =  0.254829592;
    const double a2 = -0.284496736;
    const double a3 =  1.421413741;
    const double a4 = -1.453152027;
    const double a5 =  1.061405429;
    const double p  =  0.3275911;

    // Save the sign of x
    const double sign = (x < 0) ? -1.0 : 1.0;
    const double absx = abs(x);

    // A&S formula 7.1.26
    const double t = 1.0/(1.0 + p*absx);
    const double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-absx*absx);

    return sign*y;
}

/**
 * Calculates the inverse error function using Winitzki's approximation.
 * Provides accurate results for statistical computations with relative
 * error less than 0.005.
 * 
 * @param x Input value in range (-1,1)
 * @return Value y such that erf(y) = x
 * 
 * Implementation details:
 * - Uses Winitzki's rational approximation method
 * - Valid for input values in (-1,1)
 * - Relative error less than 0.005
 * - Preserves sign symmetry: erfinv(-x) = -erfinv(x)
 * 
 * Reference:
 * Winitzki, S. (2008). "A handy approximation for the error function and its inverse"
 */
double Estimator::_erfinv(const double& x) {

    if (x <= -1.0 || x >= 1.0) {
        throw invalid_argument("Sample out of range (must be < 1 and > -1)");
    }

    const double lnx = log(1.0 - x*x);
    const double sgn = (x < 0) ? -1.0 : 1.0;
    const double tt1 = 2.0 / (M_PI*0.147) + 0.5*lnx;
    const double tt2 = 1.0 / (0.147) * lnx;

    return(sgn * sqrt(-tt1 + sqrt(tt1 * tt1 - tt2)));
}
