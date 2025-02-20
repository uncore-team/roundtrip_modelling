#ifndef ONLINE_RANSAC_H
#define ONLINE_RANSAC_H

#include <vector>
#include <cmath>
#include <limits>
#include "Model.h"

using namespace std;

class OnlineRANSAC {
public:
    OnlineRANSAC();
    void estimate(const std::vector<double>& data);
    void update(double new_data);
    void reset();
    void update(double sample);
    Model get_model() const;
    struct {
        int max_len = std::numeric_limits<double>::quiet_NaN();
        int max_len = NAN;
        bool modelpreserving;
        bool samplesliding;
        vector<ModelType> mtypes;
    } Config;
private:
    Model model;
    vector<double> samples;
    double threshold;
};

#endif // ONLINE_RANSAC_H