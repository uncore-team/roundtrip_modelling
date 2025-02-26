#include "LogNormalEstimator.h"
#include "alglib/optimization.h"

using namespace std;
using namespace alglib;

LogNormalEstimator::LogNormalEstimator() : Estimator(10) {
// Constructor
}

Model LogNormalEstimator::fit(const vector<double>& samples) {

    return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN, NAN}, {Inf, NAN}}
}

tuple<bool, GoF> LogNormalEstimator::gof(const ModelParams& params, const vector<double>& samples) { // bool previous_model
    return {true, {Inf, NaN}};
}
