/**
 * @brief Implementation of the OnlineRANSAC class for real-time outlier detection.
 * 
 * Provides methods for:
 * - Online model fitting and assessment
 * - Multiple distribution support (Log-logistic, Log-normal, Exponential)
 * - Model preservation strategies
 * - Sample sliding/forgetting mechanisms
 * - Dynamic model updating
 * 
 * Based on RANSAC (Random Sample Consensus) algorithm adapted for online processing.
 */

#include <stdexcept>
#include <tuple>
#include <optional>
#include <memory>

#include "OnlineRANSAC.h"
#include "LogLogisticEstimator.h"
#include "LogNormalEstimator.h"
#include "ExponentialEstimator.h"
#include "Model.h"

using namespace std;

/**
 * @brief Constructor for OnlineRANSAC class.
 * 
 * @param min_len Minimum number of samples needed to start modeling
 * @param max_len Maximum number of samples to maintain
 * @param model_preserving Model preservation strategy:
 *        - 0: No preservation
 *        - 1: First preserve, then rebuild
 *        - 2: First rebuild, then preserve
 * @param sample_sliding Sample window strategy:
 *        - 0: Forgetting (reset)
 *        - 1: Sliding window
 * @param data_preserving Data preservation strategy:
 *        - 0: Forgetting
 *        - 1: Preserving
 * @param model_types Vector of distribution types to try fitting
 * 
 * @throws invalid_argument if:
 *         - min_len <= 1
 *         - max_len <= min_len
 *         - model_types is empty
 */
OnlineRANSAC::OnlineRANSAC(unsigned min_len, unsigned max_len, unsigned model_preserving = 0, unsigned sample_sliding = 0, unsigned data_preserving = 0, vector<ModelType> model_types = {}) {

    // sanity checks
    if (min_len <= 1) {
        throw invalid_argument("Invalid 'min_len' param. Must be greater than 1.");
    }
    if (max_len <= min_len) {
        throw invalid_argument("Invalid 'max_len' param. Must be greater than 'min_len'.");
    }
    if (model_types.empty()) {
        throw invalid_argument("No model types param specified. Must have at least one model type.");
    }

    // save parameters
    m_min_len = min_len;
    m_max_len = max_len;
    m_model_preserving = model_preserving;
    m_sample_sliding = sample_sliding;
    m_data_preserving = data_preserving;
    m_model_types = model_types;

    // construct model estimators
    m_model_estimators = {};
    for(const ModelType& mtype: model_types) {
        switch (mtype) {
            case ModelType::LL3: m_model_estimators[mtype] = make_shared<LogLogisticEstimator>(); break;
            case ModelType::LN3: m_model_estimators[mtype] = make_shared<LogNormalEstimator>(); break;
            case ModelType::EXP: m_model_estimators[mtype] = make_shared<ExponentialEstimator>(); break;
            default:
                throw invalid_argument("Invalid model type");
        }
    }
}

/**
 * @brief Assesses data against multiple distribution models
 * 
 * Tries to fit each model type in order until finding one that fits well.
 * 
 * @param samples Vector of observations to fit
 * @return Model structure containing the first successful fit, or empty if none fit
 */
Model OnlineRANSAC::assess(const vector<double>& samples) {

    for (auto& [key, estimator] : m_model_estimators) {
        Model model = estimator->fit(samples);
        if (model.defined) {
            return model;
        }
    }
    return Model();
}

/**
 * @brief Updates the model with a new sample
 * 
 * Implements one step of the online RANSAC algorithm:
 * 1. Adds new sample to data
 * 2. Updates/preserves model based on strategy
 * 3. Maintains sample window
 * 
 * @param sample New observation to process
 * @return Exit branch indicating result:
 *         1: No model yet
 *         2: Model updated
 *         3: Model preserved
 * 
 * @throws invalid_argument if model_preserving is invalid
 */
int OnlineRANSAC::update(double sample) {

    vector<double> samples = m_state.samples;
    Model model = m_state.model;

    samples.push_back(sample); // add new sample to the list
    unsigned len = samples.size();

    // --- assessing model (algorithm core)
    int exitbranch;
    if (len < m_min_len) {
        exitbranch = 1;
    } else if (len == m_min_len) {
        model = assess(samples); // create new model (overwrites the old one)
        if(model.defined){
            exitbranch = 2;
        } else {
            exitbranch = 1;
        }
    } else {
        switch (m_model_preserving) {
            case 0: {
                model = assess(samples); // create new model (overwrites the old one)
                if (model.defined) {
                    exitbranch = 2;
                } else {
                    exitbranch = 1;
                }
                break;
            }
            case 1: {
                bool preserved = false;
                if (model.defined) {
                    Estimator::Ptr& estimator = m_model_estimators[model.type];
                    auto [reject, gof] = estimator->gof(model.params, samples); // assess existing model
                    if (!reject) {
                        model.gof = gof; // update gof if state is preserved
                        exitbranch = 3;
                        preserved = true;
                    }
                }
                if (!preserved) {
                    model = assess(samples); // create new model (overwrites the old one)
                    if (model.defined) {
                        exitbranch = 2;
                    } else {
                        exitbranch = 1;
                    }
                }
                break;
            }
            case 2: {
                model = assess(samples); // create new model
                if (model.defined) {
                    exitbranch = 2;
                } else {
                    if (model.defined) {
                        Estimator::Ptr& estimator = m_model_estimators[model.type];
                        auto [reject, gof] = estimator->gof(model.params, samples); // assess existing model
                        if (!reject) {
                            model.gof = gof; // update gof if state is updated
                            exitbranch = 3;
                        } else {
                            exitbranch = 1;
                        }
                    } else {
                        exitbranch = 1;
                    }
                }
                break;
            }
            default:
                throw invalid_argument("Invalid model-preserving");
        }
    }

    // --- preparing new state
    if (len > m_max_len) {
        samples = vector<double>(samples.end() - m_max_len, samples.end());
        // if the bounded sample violates current consensus, that will be
        // detected in the next iteration by the core algorithm
    }

    switch (exitbranch) {
        case 1: { // no model yet
            if (len > m_min_len && !m_sample_sliding) {
                samples = {sample};
            }
            else if (len >= m_min_len && !m_data_preserving) {
                samples = vector<double>(samples.end() - (m_min_len - 1), samples.end());
            }
            m_state.samples = samples;
            m_state.model = Model();
            break;
        }
        case 2: // model + consensus created now
        case 3: { // model preserved
            m_state.samples = samples;
            m_state.model = model; 
            break;
        }
        default:
            throw invalid_argument("Invalid exit branch");
    }

    return exitbranch;
}

/**
 * @brief Resets the algorithm state
 * 
 * Clears current model and samples, returning to initial state
 */
void OnlineRANSAC::reset() {
    m_state = State();
}

/**
 * @brief Returns the current model
 * 
 * @return Current model structure, or empty if no model exists
 */
Model OnlineRANSAC::get_model() {
    return m_state.model;
}
