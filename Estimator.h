#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <memory>
#include <vector>
#include <random>
#include <stdexcept>
#include <tuple>

#include "Model.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;

/**
 * @brief Abstract base class for statistical distribution fitting.
 * 
 * Provides a common interface for different distribution estimators.
 * Implements utility functions for basic statistical calculations.
 * Each derived class implements specific distribution fitting methods.
 */
class Estimator {
public:
    /** @brief Smart pointer type alias for polymorphic usage */
    using Ptr = shared_ptr<Estimator>;

    /**
     * @brief Constructor setting minimum sample size requirement
     * @param min_len Minimum number of samples needed for fitting
     */
    Estimator(unsigned min_len);

    /**
     * @brief Virtual destructor for proper cleanup of derived classes
     */
    virtual ~Estimator();

    /**
     * @brief Fits a statistical distribution to sample data
     * 
     * @param samples Vector of observations to fit
     * @return Model structure containing:
     *         - defined: true if fit succeeded
     *         - type: Distribution type identifier
     *         - params: Distribution parameters
     *         - gof: Goodness of fit statistics
     */
    virtual Model fit(const vector<double>& samples) = 0;

    /**
     * @brief Performs goodness-of-fit test on fitted distribution
     * 
     * @param params Distribution parameters to test
     * @param samples Vector of observations to test against
     * @return tuple<bool, GoF>:
     *         - bool: true if fit should be rejected
     *         - GoF: {test statistic, critical value}
     */
    virtual tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) = 0;

    /**
     * @brief Calculates the cumulative distribution function (CDF) for a single value
     * 
     * @param params Distribution parameters specific to each derived class
     * @param sample Single input value to evaluate
     * @return CDF value at the given sample point
     */
    virtual double cdf(const ModelParams& params, const double& sample) = 0;

    /**
     * @brief Calculates the CDF for a vector of values
     * 
     * @param params Distribution parameters specific to each derived class
     * @param samples Vector of input values to evaluate
     * @return Vector of CDF values corresponding to each input sample
     */
    virtual vector<double> cdf(const ModelParams& params, const vector<double>& samples) = 0;

    /**
     * @brief Calculates the probability density function (PDF) for a single value
     * 
     * @param params Distribution parameters specific to each derived class
     * @param sample Single input value to evaluate
     * @return PDF value at the given sample point
     */
    virtual double pdf(const ModelParams& params, const double& sample) = 0;

    /**
     * @brief Calculates the PDF for a vector of values
     * 
     * @param params Distribution parameters specific to each derived class
     * @param samples Vector of input values to evaluate
     * @return Vector of PDF values corresponding to each input sample
     */
    virtual vector<double> pdf(const ModelParams& params, const vector<double>& samples) = 0;

    /**
     * @brief Generates a single random value from the distribution
     * 
     * @param params Distribution parameters specific to each derived class
     * @return Single random value following the distribution
     */
    virtual double rnd(const ModelParams& params) = 0;

    /**
     * @brief Generates multiple random values from the distribution
     * 
     * @param params Distribution parameters specific to each derived class
     * @param length Number of random values to generate
     * @return Vector of random values following the distribution
     */
    virtual vector<double> rnd(const ModelParams& params, const unsigned& length) = 0;

    /**
     * @brief Calculates the expected value (mean) of the distribution
     * 
     * @param params Distribution parameters specific to each derived class
     * @return Expected value E[X] of the distribution
     */
    virtual double expectation(const ModelParams& params) = 0;

    /**
     * @brief Calculates the variance of the distribution
     * 
     * @param params Distribution parameters specific to each derived class
     * @return Variance Var[X] of the distribution
     */
    virtual double variance(const ModelParams& params) = 0;

    /**
     * @brief Calculates the mode (most frequent value) of the distribution
     * 
     * @param params Distribution parameters specific to each derived class
     * @return Mode of the distribution (highest point of PDF)
     */
    virtual double mode(const ModelParams& params) = 0;

protected:
    unsigned m_min_len;  ///< Minimum required sample size

    /**
     * @brief Random number generation members
     */
    random_device m_rnd_dev;
    mt19937 m_rnd_gen;
    uniform_real_distribution<double> m_unif_dist;  ///< Uniform distribution between 0 and 1

    /**
     * @brief Calculates the error function for a given value
     * 
     * @param x Input value
     * @return Value of the error function erf(x)
     */
    double _erf(const double& x);

    /**
     * @brief Calculates the inverse error function using numerical approximation
     * 
     * @param x Input value in range (-1,1)
     * @return Value y such that erf(y) = x
     */
    double _erfinv(const double& x);

    /**
     * @brief Utility functions for statistical calculations
     */
    double _min(const double& a, const double& b);  ///< Returns minimum of two values
    double _min(const vector<double>& samples);     ///< Returns minimum of vector
    double _max(const double& a, const double& b);  ///< Returns maximum of two values
    double _max(const vector<double>& samples);     ///< Returns maximum of vector
    double _mean(const vector<double>& samples);    ///< Calculates arithmetic mean of vector
    double _median(const vector<double>& samples);  ///< Calculates median of vector
};

#endif // ESTIMATOR_H