#include <algorithm>
#include <numeric>
#include <random>

#include <gtest/gtest.h>
#include "ExponentialEstimator.h"

using namespace std;

class ExponentialTest : public ::testing::Test {
protected:
    ExponentialEstimator estimator;
    const size_t sample_size = 5000;
    const double alpha_true = 1.0;  // location parameter
    const double beta_true = 2.0;   // rate parameter
    vector<double> samples;

    void SetUp() override {
        mt19937 gen(42); // Fixed seed for reproducibility
        exponential_distribution<double> dist(beta_true);

        samples.resize(sample_size);
        for(size_t i = 0; i < sample_size; ++i) {
            samples[i] = alpha_true + dist(gen);
        }
    }
};

TEST_F(ExponentialTest, ValidateParameterEstimation) {
    auto model = estimator.fit(samples);
    
    ASSERT_TRUE(model.defined);
    EXPECT_NEAR(model.params.alpha, alpha_true, 0.1);
    EXPECT_NEAR(model.params.beta, beta_true, 0.1);
}

TEST_F(ExponentialTest, ValidateGoodnessOfFit) {

    ModelParams params;
    params.alpha = alpha_true;
    params.beta = beta_true;

    // auto model = estimator.fit(samples);
    // auto [reject, gof] = estimator.gof(model.params, samples);
    auto [reject, gof] = estimator.gof(params, samples);

    // For properly generated data, we expect not to reject the null hypothesis
    EXPECT_FALSE(reject);
    EXPECT_LT(gof.stat, gof.thresh); // Common significance level
}

/*
TEST_F(ExponentialTest, ValidateOutlierDetection) {
    // Generate some outliers
    vector<double> data = samples;
    data.push_back(alpha_true + 10.0/beta_true); // Add obvious outlier
    
    auto model = estimator.fit(data);
    auto outliers = estimator.findOutliers(data, model.params);
    
    EXPECT_FALSE(outliers.empty());
    EXPECT_TRUE(outliers.back()); // Last point should be marked as outlier
}
*/

TEST_F(ExponentialTest, ThrowsOnInsufficientSamples) {
    vector<double> small_sample = {1.0, 2.0};
    EXPECT_THROW(estimator.fit(small_sample), invalid_argument);
}

TEST_F(ExponentialTest, ValidatePDF) {
    ModelParams params;
    params.alpha = alpha_true;
    params.beta = beta_true;
    
    // Test at multiple points
    vector<double> test_points = {1.0, 1.5, 2.0, 2.5, 3.0};
    for(double x : test_points) {
        const double z = x - alpha_true;
        if (z < 0) {
            EXPECT_EQ(estimator.pdf(params, x), 0.0);
        }
        else {
            const double expected_pdf = beta_true*exp(-beta_true*z);
            EXPECT_NEAR(estimator.pdf(params, x), expected_pdf, 1e-10);
        }
    }
}

TEST_F(ExponentialTest, ValidateCDF) {
    ModelParams params;
    params.alpha = alpha_true;
    params.beta = beta_true;

    // Test at multiple points
    vector<double> test_points = {1.0, 1.5, 2.0, 2.5, 3.0};
    for(double x : test_points) {
        const double z = x - alpha_true;
        if (z <= 0) {
            EXPECT_EQ(estimator.cdf(params, x), 0.0);
        }
        else {
            const double expected_cdf = 1.0 - exp(-beta_true*z);
            EXPECT_NEAR(estimator.cdf(params, x), expected_cdf, 1e-10);
        }
    }
}

TEST_F(ExponentialTest, ValidateVectorizedOperations) {

    ModelParams params;
    params.alpha = alpha_true;
    params.beta = beta_true;
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

TEST_F(ExponentialTest, ValidateRandomGeneration) {
    ModelParams params;
    params.alpha = alpha_true;
    params.beta = beta_true;
    const size_t n_samples = 10000;

    const vector<double> random_samples = estimator.rnd(params, n_samples);

    // Check vector size
    ASSERT_EQ(random_samples.size(), n_samples);

    // Check all samples are greater than or equal to alpha
    EXPECT_TRUE(all_of(random_samples.begin(), random_samples.end(),
                        [this](double x) { return x >= alpha_true; }));
    double accum = 0.0;
    const double n_samplesf = static_cast<double>(n_samples);

    // Statistical tests
    accum = 0.0; // Sample mean
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(n_samples > OMP_THRESH)
    #endif
    for (const double x: random_samples) { accum += x; }
    const double sample_mean = accum/n_samplesf;
    const double theoretical_mean = alpha_true + 1.0/beta_true;

    accum = 0.0; // Sample variance
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(n_samples > OMP_THRESH)
    #endif
    for (const double x: random_samples) { accum += pow(x - sample_mean, 2.0); }
    const double sample_variance = accum/(n_samplesf - 1.0);
    const double theoretical_variance = 1.0/(beta_true*beta_true);

    // Check sample statistics (with reasonable tolerances for random data)
    EXPECT_NEAR(sample_mean, theoretical_mean, 0.1);
    EXPECT_NEAR(sample_variance, theoretical_variance, 0.1);
}

TEST_F(ExponentialTest, ValidateStatistics) {

    ModelParams params;
    params.alpha = alpha_true;
    params.beta = beta_true;

    // Test expectation
    EXPECT_NEAR(estimator.expectation(params), alpha_true + 1.0/beta_true, 1e-10);

    // Test variance
    EXPECT_NEAR(estimator.variance(params), 1.0/(beta_true*beta_true), 1e-10);
}

TEST_F(ExponentialTest, ValidateBoundaryConditions) {

    // Test with edge cases
    ModelParams params;
    // params.alpha = 0.5;
    // params.beta = 1.0;
    params.alpha = alpha_true;
    params.beta = beta_true;

    EXPECT_NEAR(estimator.pdf(params, params.alpha), params.beta, 1e-10);
    EXPECT_NEAR(estimator.pdf(params, Inf), 0.0, 1e-10);

    EXPECT_NEAR(estimator.cdf(params, params.alpha), 0.0, 1e-10);
    EXPECT_NEAR(estimator.cdf(params, Inf), 1.0, 1e-10);
}
