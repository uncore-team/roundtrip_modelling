// filepath: /ransaconline2/ransaconline2/src/estimators/ExponentialEstimator.cpp
#include "ExponentialEstimator.h"
#include <stdexcept>
#include <cmath>

using namespace std;

// ExponentialEstimator class implementation
ExponentialEstimator::ExponentialEstimator() : min_len(10) {
    // Constructor
}

bool ExponentialEstimator::fit(const std::vector<double>& data, ModelParams& params) {
    // According to D'Agostino, p. 141 but corrected for the final beta to be
    // 1/mean (not correct in the book, where they say they estimate beta when 
    // actually they are estimating the mean).
    // Shifted exponential pdf in D'Agostino p. 133
    // 1/beta == mean of the distribution
    // alpha == location of the distribution
    // returns: alpha, beta, and a boolean indicating success (true) or failure (false)
    // of the fitting process 
    int len = data.size();
    if (len < min_len) {
        throw invalid_argument("Cannot fit anything with less than 10 values");
    }

    double min = *min_element(x.begin(), x.end());
    double mean = accumulate(x.begin(), x.end(), 0.0) / len;
    double mu = len * (mean - min) / (len - 1); // estimate of the (non-shifted) mean

    params.EXP.alpha = min - mu / len;
    params.EXP.beta = 1 / mu; // beta is the reciprocal of the mean (in the book they use beta as the mean)

    return true;
}

std::tuple<bool, double, double> ExponentialEstimator::assess(const vector<double>& data, Model& model) {
    // Assess the goodness of fit for the exponential model
    if (data.empty()) {
        throw std::invalid_argument("Data vector is empty.");
    }

    double statistic = 0.0;
    for (double value : data) {
        statistic += std::pow(value - beta, 2) / beta; // Chi-squared statistic
    }
    double threshold = 0.05; // Significance level

    bool reject = (statistic > threshold);
    return {reject, statistic, threshold};
}


tuple<bool, double, double> ExponentialGof(const vector<double>& x, double alpha, double beta) {
    // Based on D'Agostino p. 141: both parameters unknown. Same corrections as
    // explained in ExponentialFit() apply here.

    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distribution.");
    }

    int n = x.size();
    vector<double> xsorted = x;
    sort(xsorted.begin(), xsorted.end());
    double mu = 1 / beta; // they use beta in the book when actually they want to use the mean 
    vector<double> wsorted(n);
    transform(xsorted.begin(), xsorted.end(), wsorted.begin(), [alpha, mu](double xi) { return (xi - alpha) / mu; });
    vector<double> z(n);
    transform(wsorted.begin(), wsorted.end(), z.begin(), [](double wi) { return 1 - exp(-wi); });
    vector<double> Z = z;
    sort(Z.begin(), Z.end());

    // calculate statistic: A2 for case 3 (both parameters were deduced from
    // the same sample). This statistic measures the squared distance
    // between experimental and theoretical Zs, and, indirectly, between 
    // theoretical and experimental Xs (p. 100)
    double sumatoria = 0.0;
    for (int i = 0; i < n; ++i) {
        sumatoria += (2 * (i + 1) - 1) * log(Z[i]) + (2 * n + 1 - 2 * (i + 1)) * log(1 - Z[i]);
    }
    double A2 = -n - (1.0 / n) * sumatoria;
    // do the following since parameters come from sample (D'Agostino table 4.14)
    A2 *= (1 + 5.4 / n - 11.0 / (n * n));

    double stat = A2; // this statistic follows certain right-tailed distribution. We can set in
                      // that distribution a threshold value (in its support)
                      // corresponding to a given significance level. 
                      // Then, if the value calculated for the statistic falls 
                      // above that threshold, the hypothesis should be rejected
                      // (this is easier as the significance level grows).
                      // The p-value is the probability of the statistic distribution to
                      // produce a value of the statistic equal to or greater than the
                      // calculated one. The p-value will shrink as more strongly rejected
                      // is the null hypothesis. We do not calculate it here
                      // because the distribution of the statistic is not found in
                      // the book.
    
    // test the hypothesis 
    double thresh = 1.321; // D'Agostino table 4.14
    bool reject = (stat > thresh); // equivalently, the p-value is smaller than the significant level

    return {reject, stat, thresh};
}

pair<double, double> ExponentialFit(const vector<double>& x) {
    // According to D'Agostino, p. 141 but corrected for the final beta to be
    // 1/mean (not correct in the book, where they say they estimate beta when 
    // actually they are estimating the mean).
    // Shifted exponential pdf in D'Agostino p. 133
    // 1/beta == mean of the distribution
    // alpha == location of the distribution

    const int minlen = 10;

    int n = x.size();
    if (n < minlen) {
        throw invalid_argument("Cannot fit anything with less than 10 values");
    }

    double mi = *min_element(x.begin(), x.end());
    double mean_x = accumulate(x.begin(), x.end(), 0.0) / n;
    double mu = n * (mean_x - mi) / (n - 1); // estimate of the (non-shifted) mean
    double alpha = mi - mu / n;
    double beta = 1 / mu; // beta is the reciprocal of the mean (in the book they use beta as the mean)

    return {alpha, beta};
}