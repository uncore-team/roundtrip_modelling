#include <algorithm>
#include <numeric>
#include <random>

#include <gtest/gtest.h>
#include "LogNormalEstimator.h"

using namespace std;

class LogNormalTest : public ::testing::Test {
protected:
    LogNormalEstimator estimator;
    const size_t sample_size = 5000;

    const double gamma_true = 0.5;   // Location parameter (minimum possible value)
    const double mu_true = 1.0;      // Mean of the associated normal distribution (μ)
    const double sigma_true = 2.0;   // Standard deviation of the associated normal distribution (σ > 0)

    vector<double> samples;

    void SetUp() override {

        mt19937 gen(42); // Fixed seed for reproducibility
        lognormal_distribution<double> dist(mu_true, sigma_true);

        samples.resize(sample_size);
        for(size_t i = 0; i < sample_size; ++i) {
            samples[i] = gamma_true + dist(gen);
        }
    }
};

TEST_F(LogNormalTest, ValidateParameterEstimation) {

    auto model = estimator.fit(samples);

    ASSERT_TRUE(model.defined);
    EXPECT_NEAR(model.params.gamma, gamma_true, 0.05);
    EXPECT_NEAR(model.params.mu, mu_true, 0.05);
    EXPECT_NEAR(model.params.sigma, sigma_true, 0.05);
}

TEST_F(LogNormalTest, ValidateGoodnessOfFit) {

    ModelParams params;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

    // auto model = estimator.fit(samples);
    // auto [reject, gof] = estimator.gof(model.params, samples);
    auto [reject, gof] = estimator.gof(params, samples);
    
    // For properly generated data, we expect not to reject the null hypothesis
    EXPECT_FALSE(reject);
    EXPECT_LT(gof.stat, gof.thresh); // Common significance level
}

/*
TEST_F(LogNormalTest, ValidateOutlierDetection) {
    // Generate some outliers
    vector<double> data = samples;
    data.push_back(alpha_true + 10.0/beta_true); // Add obvious outlier
    
    auto model = estimator.fit(data);
    auto outliers = estimator.findOutliers(data, model.params);
    
    EXPECT_FALSE(outliers.empty());
    EXPECT_TRUE(outliers.back()); // Last point should be marked as outlier
}
*/

TEST_F(LogNormalTest, ThrowsOnInsufficientSamples) {

    vector<double> small_sample = {1.0, 2.0};
    EXPECT_THROW(estimator.fit(small_sample), invalid_argument);
}

TEST_F(LogNormalTest, ValidatePDF) {

    ModelParams params;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

    const double sqrt2pi = sqrt(2*M_PI);
    const double two_sigmaSq = 2*sigma_true*sigma_true;

    // Test at multiple points
    vector<double> test_points = {1.0, 1.5, 2.0, 2.5, 3.0};
    for(double x : test_points) {
        const double z = x - gamma_true;
        if (z < 0) {
            EXPECT_EQ(estimator.pdf(params, x), 0.0); // PDF debe ser 0 si x <= gamma_true
        }
        else {
            const double exponent = -pow(log(z) - mu_true, 2)/two_sigmaSq;
            const double expected_pdf = exp(exponent)/(z*sigma_true*sqrt2pi);
            EXPECT_NEAR(estimator.pdf(params, x), expected_pdf, 1e-10);
        }        
    }
}

TEST_F(LogNormalTest, ValidateCDF) {

    ModelParams params;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

    //const double sqrt2pi = sqrt(2.0*M_PI);
    const double sigma_sqrt2 = sigma_true*sqrt(2.0);

    // Test at multiple points
    vector<double> test_points = {1.0, 1.5, 2.0, 2.5, 3.0};
    for(double x : test_points) {
        const double z = x - gamma_true;
        if (z <= 0) {
            EXPECT_EQ(estimator.cdf(params, x), 0.0);
        }
        else {
            const double expected_cdf = 0.5*(1.0 + erf((log(z) - mu_true)/sigma_sqrt2));
            EXPECT_NEAR(estimator.cdf(params, x), expected_cdf, 1e-10);
        }
    }
}

TEST_F(LogNormalTest, ValidateVectorizedOperations) {

    ModelParams params;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

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

TEST_F(LogNormalTest, ValidateRandomGeneration) {

    ModelParams params;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

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
    const double theoretical_mean = gamma_true + exp(mu_true + 0.5*sigma_true*sigma_true);

    accum = 0.0; // Sample variance
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:accum) if(n_samples > OMP_THRESH)
    #endif
    for (const double x: random_samples) { accum += pow(x - sample_mean, 2.0); }
    const double sample_variance = accum/(n_samplesf - 1);
    const double theoretical_variance = (exp(sigma_true*sigma_true) - 1)*exp(2.0*mu_true + sigma_true*sigma_true);

    // Check sample mean and variance (with reasonable tolerances for random data)
    EXPECT_NEAR(sample_mean, theoretical_mean, 0.5);
    EXPECT_NEAR(sample_variance, theoretical_variance, 0.5);
}

TEST_F(LogNormalTest, ValidateStatistics) {

    ModelParams params;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

    // Test expectation
    const double theoretical_mean = gamma_true + exp(mu_true + 0.5*sigma_true*sigma_true);
    EXPECT_NEAR(estimator.expectation(params), theoretical_mean, 1e-10);

    // Test variance
    const double theoretical_variance = (exp(sigma_true*sigma_true) - 1)*exp(2.0*mu_true + sigma_true*sigma_true);
    EXPECT_NEAR(estimator.variance(params), theoretical_variance, 1e-10);
}

TEST_F(LogNormalTest, ValidateBoundaryConditions) {

    ModelParams params;
    // params.gamma = 1;
    // params.mu = 2;
    // params.sigma = 4;
    params.gamma = gamma_true;
    params.mu = mu_true;
    params.sigma = sigma_true;

    // Test with edge cases
    EXPECT_NEAR(estimator.pdf(params, params.gamma), 0.0, 1e-10);
    EXPECT_NEAR(estimator.pdf(params, Inf), 0.0, 1e-10);

    EXPECT_NEAR(estimator.cdf(params, 0.0), 0.0, 1e-10);
    EXPECT_NEAR(estimator.cdf(params, Inf), 1.0, 1e-10);
}
