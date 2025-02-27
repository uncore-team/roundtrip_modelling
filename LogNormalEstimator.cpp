#ifdef _OPENMP
    #include <omp.h>
    #pragma message("Compiling LogNormalEstimator with OpenMP support.")
#else
    #pragma message("Compiling LogNormalEstimator without OpenMP support.")
#endif
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

double LogNormalEstimator::cdf(const ModelParams& params, const double& sample) {

    double cdf;
    // double cdf = ... ;
    return cdf;
}

vector<double> LogNormalEstimator::cdf(const ModelParams& params, const vector<double>& samples) {

    // Calculate the CDF of the lognormal distribution
    unsigned len = samples.size();
    vector<double> cdf(len);
    // for (unsigned i = 0; i < len; ++i) {
    //     cdf[i] = 1.0 - exp(-beta * (samples[i] - alpha));
    // }
    return cdf;
}

double LogNormalEstimator::pdf(const ModelParams& params, const double& sample) {

    double pdf;
    // double pdf = ... ;
    return pdf;
}

vector<double> LogNormalEstimator::pdf(const ModelParams& params, const vector<double>& samples) {

    // Calculate the PDF of the lognormal distribution
    unsigned len = samples.size();
    vector<double> pdf(len);
    // code for calculating the pdf here
    // for (unsigned i = 0; i < len; ++i) {
    //     pdf[i] = ...;
    // } 
    return pdf;
}
    
double LogNormalEstimator::rnd(const ModelParams& params) {

    double p = m_unif_dist(m_rnd_gen);
    // code for generating a ramdom number here
    double rnd = p;
    //
    return rnd;
}    

vector<double> LogNormalEstimator::rnd(const ModelParams& params, const unsigned& length) {

    vector<double> samples(length);
    // for (unsigned i = 0; i < length; ++i) {
    //     double p = m_unif_dist(m_rnd_gen);
    //     samples[i] = -log(1.0 - p) / beta + alpha;
    // }
    return samples;
}    
