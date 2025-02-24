#include "ExponentialEstimator.h"

using namespace std;

// ExponentialEstimator class implementation

ExponentialEstimator::ExponentialEstimator() : m_min_len(10) {
// Constructor
}

Model ExponentialEstimator::fit(const vector<double>& samples) {
// According to D'Agostino, p. 141 but corrected for the final beta to be
// 1/mean (not correct in the book, where they say they estimate beta when 
// actually they are estimating the mean).
// Shifted exponential pdf in D'Agostino p. 133
// 1/beta == mean of the distribution
// alpha == location of the distribution
// returns: alpha, beta, and a boolean indicating success (true) or failure (false)
// of the fitting process

    int len = samples.size();

    // sanity check
    if (len < m_min_len) {
        throw invalid_argument("Cannot fit anything with less than 10 values");
    }

    double min = *min_element(samples.begin(), samples.end());
    double mean = accumulate(samples.begin(), samples.end(), 0.0) / len;
    double mu = len * (mean - min) / (len - 1); // estimate of the (non-shifted) mean

    ModelParams params;
    params.alpha = min - mu / len;
    params.beta = 1 / mu; // beta is the reciprocal of the mean (in the book they use beta as the mean)

    auto [reject, gof] = this->gof(params, samples);

    if (!reject) {
        return {true, ModelType::EXP, params, gof};
    }
    else { 
        return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN}, {NAN, NAN}}
    }
}

tuple<bool, GoF> ExponentialEstimator::gof(const ModelParams& params, const vector<double>& samples) {
// Based on D'Agostino p. 141: both parameters unknown. Same corrections as
// explained in fit() apply here.

    double mu = 1 / params.beta; // they use beta in the book when actually they want to use the mean
    double alpha = params.alpha;
    double beta = params.beta;
    double thresh = 1.321; // D'Agostino table 4.14
    int len = samples.size();

    // sanity check
    if (beta <= 0) {
        throw invalid_argument("Invalid beta for exponential distribution.");
    }
    if (len < m_min_len) {
        throw invalid_argument("Number of samples is not enough or is zero.");
    }

    // prepare data for the test
    vector<double> data(len);
    transform(
        samples.begin(), samples.end(),
        data.begin(),
        [alpha, mu](double sample) { return 1 - exp(-(sample - alpha) / mu); }
    );
    sort(data.begin(), data.end());

    // calculate statistic: A2 for case 3 (both parameters were deduced from
    // the same sample). This statistic measures the squared distance
    // between experimental and theoretical Zs, and, indirectly, between 
    // theoretical and experimental Xs (p. 100)
    double accum = 0.0;
    for (int i = 0; i < len; ++i) {
        accum += (2 * (i + 1) - 1) * log(data[i]) + (2 * len + 1 - 2 * (i + 1)) * log(1 - data[i]);
    }

    double A2 = -len - (1.0 / len) * accum;

    // do the following since parameters come from sample (D'Agostino table 4.14)
    A2 *= (1 + 5.4 / len - 11.0 / (len * len));

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
    bool reject = (stat > thresh); // equivalently, the p-value is smaller than the significant level

    return {reject, {stat, thresh}};
}
