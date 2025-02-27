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
 * @brief Calculates arithmetic mean of sample data
 * @param samples Vector of observations
 * @return Arithmetic mean of the samples
 */
double Estimator::_mean(const vector<double>& samples) {
    unsigned len = samples.size();
    double mean = accumulate(samples.begin(), samples.end(), 0.0) / len;
    return mean;
}

/**
 * @brief Calculates median of sample data
 * @param samples Vector of observations (copied for sorting)
 * @return Median value:
 *         - For odd n: middle value
 *         - For even n: average of two middle values
 */
double Estimator::_median(vector<double> samples) {
    unsigned len = samples.size();
    sort(samples.begin(), samples.end());
    if (len % 2 != 0) 
        return samples[len/2];
    else
        return (samples[(len-1)/2] + samples[len/2])/2.0;
}

/**
 * @brief Returns minimum of two values
 * @param a First value to compare
 * @param b Second value to compare
 * @return Smaller of the two values
 */
double Estimator::_min(const double& a, const double& b) {
    return ((a) < (b) ? a : b);
}

/**
 * @brief Finds minimum value in sample data
 * @param samples Vector of observations
 * @return Minimum value in the vector
 */
double Estimator::_min(const vector<double>& samples) {
    double min = *min_element(samples.begin(), samples.end());
    return min;
}

/**
 * @brief Returns maximum of two values
 * @param a First value to compare
 * @param b Second value to compare
 * @return Larger of the two values
 */
double Estimator::_max(const double& a, const double& b) {
    return ((a) > (b) ? a : b);
}

/**
 * @brief Finds maximum value in sample data
 * @param samples Vector of observations
 * @return Maximum value in the vector
 */
double Estimator::_max(const vector<double>& samples) {
    double max = *max_element(samples.begin(), samples.end());
    return max;
}
