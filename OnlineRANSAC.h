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

class OnlineRANSAC {
public:
    OnlineRANSAC(unsigned min_len, unsigned max_len, unsigned model_preserving, bool sample_sliding, bool data_preserving, vector<ModelType> model_types);
    void reset();
    int update(double sample);
    Model get_model();

private:
    State m_state;
    unsigned m_min_len;
    unsigned m_max_len;
    unsigned m_model_preserving;
    bool m_sample_sliding;
    bool m_data_preserving;
    vector<ModelType> m_model_types;
    map<ModelType, std::shared_ptr<Estimator>> m_model_estimators;

    Model assess(const vector<double>& samples);
};

#endif // ONLINE_RANSAC_H