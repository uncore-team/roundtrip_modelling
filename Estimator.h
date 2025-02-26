#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <memory>
#include <vector>
#include <tuple>

#include "Model.h"

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
    virtual std::tuple<bool, GoF> gof(const ModelParams& params, const vector<double>& samples) = 0;

protected:
    unsigned m_min_len;  ///< Minimum required sample size

    /**
     * @brief Utility functions for statistical calculations
     */
    double _min(const double& a, const double& b);     ///< Returns minimum of two values
    double _min(const vector<double>& samples);        ///< Returns minimum of vector
    double _max(const double& a, const double& b);     ///< Returns maximum of two values
    double _max(const vector<double>& samples);        ///< Returns maximum of vector
    double _mean(const vector<double>& samples);       ///< Calculates arithmetic mean of vector
    double _median(vector<double> samples);            ///< Calculates median of vector
};

#endif // ESTIMATOR_H