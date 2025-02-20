#include <gtest/gtest.h>
#include "exponential_model.h"

class ExponentialModelTest : public ::testing::Test {
protected:
    ExponentialModel model;

    void SetUp() override {
        // Initialize the model if necessary
    }
};

TEST_F(ExponentialModelTest, TestFit) {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    model.fit(data);
    EXPECT_GT(model.getLambda(), 0); // Check that the estimated lambda is positive
}

TEST_F(ExponentialModelTest, TestProbability) {
    model.fit({1.0, 2.0, 3.0});
    double probability = model.probability(2.0);
    EXPECT_GT(probability, 0); // Check that the probability is positive
}

TEST_F(ExponentialModelTest, TestInvalidData) {
    EXPECT_THROW(model.fit({}), std::invalid_argument); // Check that fitting with empty data throws an exception
}