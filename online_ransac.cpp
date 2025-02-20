#include "online_ransac.h"
#include "model/exponential_model.h"
#include "utils/math_utils.h"

OnlineRANSAC::OnlineRANSAC() {
    // Initialize parameters and state

}

void OnlineRANSAC::estimate(const std::vector<double>& data) {
// Fit the exponential model to the provided data
}

//pair<int, State> AlgOnlineransac(const Params& parms, State state, double rt) {

void OnlineRANSAC::update(double sample) {
// Do one step of the on-line RANSAC algorithm for modelling round-trip times.
//
// PARMS -> a struct with the parameters of the algorithm:
//           .s -> minimum number of round-trip times to start modelling.
//           .N -> maximum number of round-trip times to be modelled. 
//                 If NaN, no maximum.
//           .mtypes -> a cell with as many strings as model types we wish
//                      to use for modelling, in the same order we wish to
//                      try them at modelling. See ModelCreate().
//           .modelpreserving -> either 0 (no preserving), 1 (first
//                               preserve, then rebuild), or 2 (first
//                               rebuild, then preserve).
//           .samplesliding -> either 0 (forgetting) or 1 (sliding).
//           .datapreserving -> either 0 (forgetting) or 1 (preserving).
// STATE -> a struct with the current state of the algorithm:
//           .sample -> current round-trip sample used for modelling.
//           .ind0 -> index within the global scenario of the first 
//                    round-trip time in the sample used
//                    so far.
//           .ind1 -> index of the last round-trip time in the sample used
//                    so far, within the complete scenario.
//           .model -> the model currently active (see ModelCreate) or empty
//                     if none. The model has been assessed for the SAMPLE.
//           .consensus -> info about the validity of the consensus set (see
//                         ModelGof) or empty if none.
//          Initially, this struct must contain an empty sample, model and
//          consensus, and .ind* == NaN.
// RT -> the current round-trip time, to be processed in this call.
//
// EXITBRANCH <- reason of exit:
//               '1' <- no model yet
//               '2' <- model updated
//               '3' <- model preserved
// NEWSTATE <- updated STATE after processing the round-trip time RT.
    
    // --- sanity checks
    if (parms.s <= 1) {
        throw invalid_argument("Invalid params.s");
    }
    if (parms.mtypes.empty()) {
        throw invalid_argument("No model specified");
    }

    // --- preparation of data
    vector<double> S = state.sample;
    S.push_back(rt);
    int ind0 = (state.ind0 < 0) ? 0 : state.ind0;
    int ind1 = (state.ind0 < 0) ? 0 : state.ind0 + state.sample.size();
    int ns = S.size();
    State newstate = state;

    // --- core of the algorithm
    int exitbranch;
    optional<Model> m;
    optional<GofInfo> g;

    if (ns < parms.s) {
        exitbranch = 1;
    } else if (ns == parms.s) {
        tie(m, g) = ModelAssess(S, ind0, ind1, parms.mtypes); // create new model
        if (m && g) {
            exitbranch = 2;
        } else {
            exitbranch = 1;
        }
    } else {
        switch (parms.modelpreserving) {
            case 0: {
                tie(m, g) = ModelAssess(S, ind0, ind1, parms.mtypes); // create new model
                if (m && g) {
                    exitbranch = 2;
                } else {
                    if (!parms.samplesliding) {
                        S = {rt};
                        ind0 = ind1;
                    }
                    exitbranch = 1;
                }
                break;
            }
            case 1: {
                bool preserved = false;
                if (state.model.defined) {
                    tie(ignore, g) = ModelAssess(S, ind0, ind1, state.model); // assess existing model
                    if (g) {
                        newstate.consensus = *g; // newstate.model is preserved
                        preserved = true;
                        exitbranch = 3;
                    }
                }
                if (!preserved) {
                    tie(m, g) = ModelAssess(S, ind0, ind1, parms.mtypes); // create new model
                    if (m && g) {
                        exitbranch = 2;
                    } else {
                        if (!parms.samplesliding) {
                            S = {rt};
                            ind0 = ind1;
                        }
                        exitbranch = 1;
                    }
                }
                break;
            }
            case 2: {
                auto oldmodel = state.model;
                tie(m, g) = ModelAssess(S, ind0, ind1, parms.mtypes); // create new model
                if (m && g) {
                    exitbranch = 2;
                } else {
                    if (oldmodel.defined) { // oldmodel == state.model
                        tie(ignore, g) = ModelAssess(S, ind0, ind1, state.model); // assess existing model
                        if (g) {
                            newstate.consensus = *g;
                            exitbranch = 3;
                        } else {
                            if (!parms.samplesliding) {
                                S = {rt};
                                ind0 = ind1;
                            }
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
    switch (exitbranch) {
        case 1: { // no model yet
            ns = S.size(); // may have varied due to sliding; ind0,ind1 updated already
            if (!parms.datapreserving && ns >= parms.s) { // parms.s > 2 here
                // preserve only the newest s-1 values
                newstate.sample = vector<double>(S.end() - (parms.s - 1), S.end());
                newstate.ind0 = ind1 - (parms.s - 1) + 1; // ind1 points to the last one in S right before this
                newstate.ind1 = ind1;
            } else {
                newstate.sample = S;
                newstate.ind0 = ind0;
                newstate.ind1 = ind1;
            }
            newstate.model = Model();
            newstate.consensus = GofInfo();
            break;
        }
        case 2: { // model + consensus created now
            newstate.sample = S;
            newstate.ind0 = ind0;
            newstate.ind1 = ind1;
            newstate.model = move(*m); // m, g calculated in the corresponding case
            newstate.consensus = move(*g);
            break;
        }
        case 3: { // model preserved
            newstate.sample = S;
            newstate.ind0 = ind0;
            newstate.ind1 = ind1;
            // model and consensus already updated in the corresponding case
            break;
        }
        default:
            throw invalid_argument("Invalid exit branch");
    }

    // --- bounding gathered sample
    ns = newstate.sample.size();
    if (parms.N >= 0 && ns > parms.N) {
        newstate.sample = vector<double>(newstate.sample.end() - parms.N, newstate.sample.end());
        newstate.ind0 = newstate.ind1 - parms.N + 1;
        // if the bounded sample violates current consensus, that will be
        // detected in the next iteration by the core algorithm
    }
    if (newstate.ind1 - newstate.ind0 + 1 != newstate.sample.size()) {
        throw invalid_argument("Invalid indexes of sample");
    }

    return {exitbranch, newstate};
}

void OnlineRANSAC::reset() {
    // Reset the model and state
    model.reset();
    current_state = State(); // Assuming State is a struct that holds the current state
}

const ExponentialModel& OnlineRANSAC::getModel() const {
    return model;
}

const State& OnlineRANSAC::getState() const {
    return current_state;
}