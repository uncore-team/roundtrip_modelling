#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <optional>
#include "estimators/LL3Estimator.h"
#include "estimators/LN3Estimator.h"
#include "estimators/BernoulliEstimator.h"
#include "estimators/ExponentialEstimator.h"
#include "types/Model.h"

using namespace std;

int main() {
    // Initialize parameters for the estimators
    Params parms = {20, -1, {"LL3"}, 0, 1, 0};

    // State initialization
    State state = {{}, -1, -1, {}, {}};

    // Read RTTs from a file
    string filename = "rtts.txt";
    vector<double> rtts = readRTTsFromFile(filename);

    // Create estimator objects
    LL3Estimator ll3Estimator;
    LN3Estimator ln3Estimator;
    BernoulliEstimator bernoulliEstimator;
    ExponentialEstimator exponentialEstimator;

    // Process the RTTs with the LL3 estimator
    auto ll3Results = ll3Estimator.performEstimation(rtts, parms, state);

    // Process the RTTs with the LN3 estimator
    auto ln3Results = ln3Estimator.performEstimation(rtts, parms, state);
    
    // Process the RTTs with the Bernoulli estimator
    auto bernoulliResults = bernoulliEstimator.performEstimation(rtts, parms, state);
    
    // Process the RTTs with the Exponential estimator
    auto exponentialResults = exponentialEstimator.performEstimation(rtts, parms, state);

    // Output the results (example for LL3)
    for (const auto& result : ll3Results) {
        cout << "LL3 Index: " << result.index << ", Performance: " << result.perf << endl;
    }

    return 0;
}