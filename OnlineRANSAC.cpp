#include <stdexcept>
#include <tuple>
#include <optional>
#include <memory>

#include "OnlineRANSAC.h"
#include "ExponentialEstimator.h"
#include "Model.h"

using namespace std;

OnlineRANSAC::OnlineRANSAC(unsigned min_len, unsigned max_len, unsigned model_preserving = 0, bool sample_sliding = false, bool data_preserving = false, vector<ModelType> model_types = {}) {
//           min_len -> minimum number of round-trip times to start modelling.
//           max_len -> maximum number of round-trip times to be modelled. 
//                      If NaN, no maximum.
//           model_preserving -> either 0 (no preserving), 1 (first
//                               preserve, then rebuild), or 2 (first
//                               rebuild, then preserve).
//           sample_sliding -> either 0 (forgetting) or 1 (sliding).
//           data_preserving -> either 0 (forgetting) or 1 (preserving).
//           model_estimators -> a map with as many ModelType as we wish
//                               to use for modelling, in the same
//                               order we wish to try them at modelling.

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
            case ModelType::LL3: /* m_model_estimators[mtype] = make_shared<LL3Estimator>();*/ break;
            case ModelType::LN3: /* m_model_estimators[mtype] = make_shared<LN3Estimator>();*/ break;
            case ModelType::EXP: m_model_estimators[mtype] = make_shared<ExponentialEstimator>(); break;
            default:
                throw invalid_argument("Invalid model type");
        }
    }
}

Model OnlineRANSAC::assess(const vector<double>& samples) {
// If MODELS is a vector with some model types, select the first one that 
// fits ok to the data in S and is assessed with GoF; if models
// is a model, just assesses it with GoF.
//
// samples -> samples to fit/gof.
// MODELS -> either a vector with as many entries as model types, each one a 
//           text, or a model already defined (see ModelCreate).
//
// M <- model fitted/gofed to the data or a copy of MODELS if
//      MODELS is a model. Empty in any case the assessment fails.

    for (auto& [key, estimator] : m_model_estimators) {
        Model model = estimator->fit(samples);
        if (model.defined) {
            return model;
        }
    }
    return Model();
}

int OnlineRANSAC::update(double sample) {
// Do one step of the on-line RANSAC algorithm for modelling round-trip times.
//
// sample -> the current round-trip time, to be processed in this call.
//
// EXITBRANCH <- reason of exit:
//               '1' <- no model yet
//               '2' <- model updated
//               '3' <- model preserved

    // --- preparation of data
    vector<double> samples = m_state.samples;

    // --- core of the algorithm
    int exitbranch;
    Model model;

    samples.push_back(sample); // add new sample to the list
    int len = samples.size();

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
                model = m_state.model;
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
                    model = m_state.model;
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
    len = samples.size();
    if (len > m_max_len) {
        samples = vector<double>(samples.end() - m_max_len, samples.end());
        // if the bounded sample violates current consensus, that will be
        // detected in the next iteration by the core algorithm
    }

    switch (exitbranch) {
        case 1: { // no model yet
            if (len >= m_min_len) {
                if (!m_sample_sliding) {
                    samples = {sample};
                }
                else if (!m_data_preserving){
                    samples = vector<double>(samples.end() - (m_min_len - 1), samples.end());
                }
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

void OnlineRANSAC::reset() {
    m_state = State();
}

Model OnlineRANSAC::get_model() {
    return m_state.model;
}
