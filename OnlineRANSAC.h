#ifndef ONLINE_RANSAC_H
#define ONLINE_RANSAC_H

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "Model.h"
#include "Estimator.h"

using namespace std;

/**
 * @brief Class for online outlier detection using RANSAC algorithm.
 * 
 * Implements a real-time adaptation of RANSAC with:
 * - Multiple distribution model support
 * - Flexible model preservation strategies
 * - Configurable sample window management
 * - Dynamic model updating and assessment
 * 
 * The class maintains state information about:
 * - Current model and its parameters
 * - Sample window and its management
 * - Model preservation strategy
 */
class OnlineRANSAC {
public:
    /**
     * @brief Constructor initializing RANSAC parameters
     * 
     * @param min_len Minimum samples needed to start modeling
     * @param max_len Maximum samples to maintain in window
     * @param model_preserving Preservation strategy (0-2)
     * @param sample_sliding Window strategy (0-1)
     * @param data_preserving Data retention strategy (0-1)
     * @param model_types Vector of distributions to try fitting
     */
    OnlineRANSAC(unsigned min_len, unsigned max_len, 
                 unsigned model_preserving, unsigned sample_sliding, 
                 unsigned data_preserving, vector<ModelType> model_types);

    /**
     * @brief Resets algorithm state to initial conditions
     */
    void reset();

    /**
     * @brief Processes a new sample and updates model
     * @param sample New observation to process
     * @return Update result code (1-3)
     */
    int update(double sample);

    /**
     * @brief Returns current model state
     * @return Model structure or empty if no model exists
     */
    Model get_model();

private:
    State m_state;                    ///< Current algorithm state
    unsigned m_min_len;               ///< Minimum required samples
    unsigned m_max_len;               ///< Maximum window size
    unsigned m_model_preserving;      ///< Preservation strategy
    unsigned m_sample_sliding;        ///< Window strategy
    unsigned m_data_preserving;       ///< Data retention strategy
    vector<ModelType> m_model_types;  ///< Supported distributions
    map<ModelType, std::shared_ptr<Estimator>> m_model_estimators;  ///< Distribution estimators

    /**
     * @brief Attempts to fit models to sample data
     * @param samples Vector of observations to fit
     * @return Best fitting model or empty if none fit
     */
    Model assess(const vector<double>& samples);
};

#endif // ONLINE_RANSAC_H