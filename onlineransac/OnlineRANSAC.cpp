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

#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>

#include "OnlineRANSAC.h"
#include "Model.h"

using namespace std;

/**
 * Constructor for OnlineRANSAC class.
 * Initializes the online outlier detection algorithm with specified parameters.
 * 
 * @param min_len Minimum number of samples needed to start modeling (> 1)
 * @param max_len Maximum number of samples to maintain (> min_len)
 * @param model_preserving Model preservation strategy:
 *        - 0: No preservation (always rebuild)
 *        - 1: First preserve, then rebuild if preservation fails
 *        - 2: First rebuild, then preserve if rebuild fails
 * @param sample_sliding Sample window strategy:
 *        - 0: Reset window on model failure
 *        - 1: Maintain sliding window
 * @param data_preserving Data preservation strategy:
 *        - 0: Reset to minimum length on failure
 *        - 1: Preserve all valid samples
 * @param model_types Vector of distribution types to try fitting
 * 
 * Implementation details:
 * - Validates all input parameters
 * - Initializes distribution estimators based on model types
 * - Sets up initial empty state
 * 
 * @throws invalid_argument if:
 *         - min_len <= 1
 *         - max_len <= min_len
 *         - model_types is empty
 *         - Invalid model type specified
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
    for(const ModelType& mtype: model_types) {
        m_model_estimators[mtype] = Estimator::create(mtype);
    }
}

/**
 * Assesses data against multiple distribution models.
 * Attempts to fit each configured distribution type in sequence.
 * 
 * @param samples Vector of observations to fit
 * @return Model structure containing:
 *         - First successful fit if any model fits
 *         - Empty model if no distribution fits well
 * 
 * Implementation details:
 * - Tries each model type in order specified at construction
 * - Uses Anderson-Darling test for goodness of fit
 * - Returns first model that passes goodness of fit test
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
 * Updates the model with a new sample.
 * Core implementation of the online RANSAC algorithm.
 * 
 * @param sample New observation to process
 * @return Exit branch indicating result:
 *         1: No valid model (insufficient data or no fit)
 *         2: New model created and validated
 *         3: Existing model preserved
 * 
 * Implementation details:
 * - Adds new sample to data window
 * - Handles initialization phase (len < min_len)
 * - Implements model preservation strategy
 * - Maintains sample window size
 * - Updates algorithm state
 * 
 * @throws invalid_argument if:
 *         - model_preserving value is invalid
 *         - invalid exit branch encountered
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
                    auto [reject, gof] = estimator->gof(model.params, samples, true); // assess existing model
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
                        auto [reject, gof] = estimator->gof(model.params, samples, true); // assess existing model
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
 * Resets the algorithm state to initial conditions.
 * Clears all current data and models.
 * 
 * Implementation details:
 * - Clears sample window
 * - Resets current model
 * - Maintains configuration parameters
 */
void OnlineRANSAC::reset() {
    m_state = State();
}

/**
 * Returns the current model state.
 * 
 * @return Current model structure containing:
 *         - Distribution type and parameters if model exists
 *         - Empty model if no valid model exists
 * 
 * Implementation details:
 * - Returns copy of internal model state
 * - Does not modify algorithm state
 */
Model OnlineRANSAC::get_model() {
    return m_state.model;
}

/**
 * Prints the current fitted model parameters.
 * Displays model type and parameters in human-readable format.
 * 
 * Implementation details:
 * - Uses fixed precision (6 decimal places)
 * - Shows distribution-specific parameters
 * - Handles undefined models
 * - Supports all configured model types
 * 
 * @throws out_of_range if unknown model type encountered
 */   
void OnlineRANSAC::print_model() {

    const Model& model = m_state.model;

    if (!model.defined) {
        cout << "\nModel Type: model NOT defined." << endl;
        return;
    }

    cout << fixed << setprecision(6);  // Configurar formato para números flotantes
    cout << "\nModel Type: ";
    string type;
    switch(model.type) {
        case ModelType::LL3:
            cout << "LL3\n\ta: " << model.params.a
            << "\n\tb: " << model.params.b
            << "\n\tc: " << model.params.c;
            break;
        case ModelType::LN3:
            cout << "LN3\n\tgamma: " << model.params.gamma
                << "\n\tmu: " << model.params.mu
                << "\n\tsigma: " << model.params.sigma;
            break;
        case ModelType::EXP:
            cout << "EXP\n\talpha: " << model.params.alpha
                 << "\n\tbeta: " << model.params.beta;
            break;
        case ModelType::None:
            cout << "Model is not defined.";
            break;
        default:
            throw out_of_range("Unknown model type.");
    }
    cout << endl;
}
