#ifndef MODEL_H
#define MODEL_H

#include <limits>
#include <cmath>

enum ModelType {
    None = 0,
    LL3 = 1,
    LN3 = 2,
    EXP = 3
};

// Structure to hold model coefficients and type
union ModelParams {
    struct LL3 { // LogLogistic (3-params) model
        double a = NAN;
        double b = NAN;
        double c = NAN;
    };
    struct LN3 { //  LogNormal (3-params) model
        double gamma = NAN;  // Gamma
        double mu = NAN;     // Mean
        double sigma = NAN;  // Standard deviation
    };
    struct EXP { // Exponential model
        double alpha = NAN;
        double beta = NAN;
    };
};

struct GoF { // Goodness of fit
    double stat;  // Statistic
    double thresh;  // Threshold
};

// Structure to represent a statistical model
struct Model {
    bool defined = false;  // Flag to indicate if the model is defined
    ModelType type = ModelType::None;  // Type of the model (e.g., "LL3", "LN3", "EXP")
    ModelParams params;  // Coefficients of the model
    GoF gof;  // Goodness of fit
};

struct State {
    vector<double> sample;
    Model model;
};

#endif // MODEL_H