#include <algorithm>
#include <numeric>
#include <random>

#include <gtest/gtest.h>
#include "LogLogisticEstimator.h"

using namespace std;

class LogLogisticTest : public ::testing::Test {
protected:
    LogLogisticEstimator estimator;
    const size_t sample_size = 5000;

    const double gamma_true = 0.5;   // Location parameter (minimum possible value)
    const double alpha_true = 2.0;   // scale parameter and is also the median of the distribution
    const double beta_true = 3.0;    // 1 / c shape

    vector<double> samples;

    void SetUp() override {

        mt19937 gen(42); // Fixed seed for reproducibility
        uniform_real_distribution<double> dist(0.0, 1.0);

        samples.resize(sample_size);
        for(size_t i = 0; i < sample_size; ++i) {
            const double u = dist(gen);
            samples[i] = gamma_true + alpha_true*pow(u/(1.0 - u), 1.0/beta_true);
        }
    
        // write_data("ll3_samples.txt", samples);
    }
};

TEST_F(LogLogisticTest, ValidateParameterEstimation) {

    auto model = estimator.fit(samples);

    ASSERT_TRUE(model.defined);
    EXPECT_NEAR(model.params.a, gamma_true, 0.05);
    EXPECT_NEAR(model.params.b, alpha_true, 0.05);
    EXPECT_NEAR(model.params.c, 1.0/beta_true, 0.05);
}

TEST_F(LogLogisticTest, ValidateGoodnessOfFit) {

    ModelParams params;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    // auto model = estimator.fit(samples);
    // auto [reject, gof] = estimator.gof(model.params, samples);
    auto [reject, gof] = estimator.gof(params, samples);

    // For properly generated data, we expect not to reject the null hypothesis
    EXPECT_FALSE(reject);
    EXPECT_LT(gof.stat, gof.thresh); // Common significance level
}

/*
TEST_F(LogLogisticTest, ValidateOutlierDetection) {
    // Generate some outliers
    vector<double> data = samples;
    data.push_back(alpha_true + 10.0/beta_true); // Add obvious outlier
    
    auto model = estimator.fit(data);
    auto outliers = estimator.findOutliers(data, model.params);
    
    EXPECT_FALSE(outliers.empty());
    EXPECT_TRUE(outliers.back()); // Last point should be marked as outlier
}
*/

TEST_F(LogLogisticTest, ThrowsOnInsufficientSamples) {

    vector<double> small_sample = {1.0, 2.0};
    EXPECT_THROW(estimator.fit(small_sample), invalid_argument);
}

TEST_F(LogLogisticTest, ValidatePDF) {

    ModelParams params;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    // Test at multiple points
    vector<double> test_points = {1.0, 1.5, 2.0, 2.5, 3.0};
    for(double x : test_points) {
        const double z = x - gamma_true;
        if (z < 0) {
            EXPECT_EQ(estimator.pdf(params, x), 0.0); // PDF should be > 0 if x >= gamma_true
        }
        else {
            const double base = (x - gamma_true)/alpha_true;
            const double expected_pdf = (beta_true/alpha_true)*pow(base, (beta_true - 1.0))/
                                            pow((1.0 + pow(base, beta_true)), 2.0);
            EXPECT_NEAR(estimator.pdf(params, x), expected_pdf, 1e-10);
        }        
    }
}

TEST_F(LogLogisticTest, ValidateCDF) {

    ModelParams params;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    // Test at multiple points
    vector<double> test_points = {1.0, 1.5, 2.0, 2.5, 3.0};
    for(double x : test_points) {
        const double z = x - gamma_true;
        if (z <= 0) {
            EXPECT_EQ(estimator.cdf(params, x), 0.0);
        }
        else {
            const double expected_cdf = 1.0/(1.0 + pow(alpha_true/(x - gamma_true), beta_true));
            EXPECT_NEAR(estimator.cdf(params, x), expected_cdf, 1e-10);
        }
    }
}

TEST_F(LogLogisticTest, ValidateVectorizedOperations) {

    ModelParams params;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    vector<double> test_points = {1.5, 2.0, 2.5};
    
    auto cdfs = estimator.cdf(params, test_points);
    auto pdfs = estimator.pdf(params, test_points);
    
    ASSERT_EQ(cdfs.size(), test_points.size());
    ASSERT_EQ(pdfs.size(), test_points.size());
    
    for(size_t i = 0; i < test_points.size(); ++i) {
        EXPECT_NEAR(cdfs[i], estimator.cdf(params, test_points[i]), 1e-10);
        EXPECT_NEAR(pdfs[i], estimator.pdf(params, test_points[i]), 1e-10);
    }
}

TEST_F(LogLogisticTest, ValidateRandomGeneration) {

    ModelParams params;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    const size_t n_samples = 10000;
    const vector<double> random_samples = estimator.rnd(params, n_samples);

    // Check vector size
    ASSERT_EQ(random_samples.size(), n_samples);

    // Check all samples are greater than or equal to gamma
    EXPECT_TRUE(all_of(random_samples.begin(), random_samples.end(),
                        [this](double x) { return x >= gamma_true; }));

    double accum = 0.0;
    const double n_samplesf = static_cast<double>(n_samples);

    accum = 0.0; // Sample mean
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(n_samples > OMP_THRESH)
    #endif
    for (const double x: random_samples) { accum += x; }
    const double sample_mean = accum/n_samplesf;
    const double theoretical_mean = gamma_true + (alpha_true*M_PI)/(beta_true*sin(M_PI/beta_true));

    accum = 0.0; // Sample variance
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(n_samples > OMP_THRESH)
    #endif
    for (const double x: random_samples) { accum += pow(x - sample_mean, 2.0); }
    const double sample_variance = accum/(n_samplesf - 1);
    const double b = (M_PI/beta_true);
    const double theoretical_variance = 
            alpha_true*alpha_true*(2.0*b/sin(2.0*b) - pow((b/sin(b)), 2.0));

    // Check sample mean and variance (with reasonable tolerances for random data)
    EXPECT_NEAR(sample_mean, theoretical_mean, 0.5);
    EXPECT_NEAR(sample_variance, theoretical_variance, 0.5);
}

TEST_F(LogLogisticTest, ValidateStatistics) {

    ModelParams params;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    // Test expectation
    if (beta_true <= 1.0)
        EXPECT_EQ(estimator.expectation(params), NaN);
    else {
        const double theoretical_mean = gamma_true + 
                                            (alpha_true*M_PI)/(beta_true*sin(M_PI/beta_true));
        EXPECT_NEAR(estimator.expectation(params), theoretical_mean, 1e-10);
    }

    // Test variance
    if (beta_true <= 2.0)
        EXPECT_EQ(estimator.variance(params), NaN);
    else
    {
        const double b = (M_PI/beta_true);
        const double theoretical_variance = 
                alpha_true*alpha_true*(2.0*b/sin(2.0*b) - pow((b/sin(b)), 2.0));
        EXPECT_NEAR(estimator.variance(params), theoretical_variance, 1e-10);
    }

    // Test mode
    if (beta_true <= 1.0)
        EXPECT_EQ(estimator.mode(params), NaN);
    else {
        const double theoretical_mode = gamma_true + alpha_true*pow((beta_true - 1)/(beta_true + 1), 1.0/beta_true);
        EXPECT_NEAR(estimator.mode(params), theoretical_mode, 1e-10);
    }
}

TEST_F(LogLogisticTest, ValidateBoundaryConditions) {

    ModelParams params;
    // params.gamma = 1;
    // params.mu = 2;
    // params.sigma = 4;
    params.a = gamma_true;
    params.b = alpha_true;
    params.c = 1/beta_true;

    // Test with edge cases
    EXPECT_NEAR(estimator.pdf(params, gamma_true), 0.0, 1e-10);
    EXPECT_NEAR(estimator.pdf(params, Inf), 0.0, 1e-10);

    EXPECT_NEAR(estimator.cdf(params, 0.0), 0.0, 1e-10);
    EXPECT_NEAR(estimator.cdf(params, Inf), 1.0, 1e-10);
}
