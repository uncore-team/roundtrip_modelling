#ifndef MODEL_H
#define MODEL_H

#include <limits>

using namespace std;

const double eps = numeric_limits<double>::epsilon();
const double Inf = numeric_limits<double>::max();
const double NaN = numeric_limits<double>::quiet_NaN();

#define MIN(a,b) ((a) < (b) ? a : b) 
#define MAX(a,b) ((a) > (b) ? a : b)

enum ModelType {
    None = 0,
    LL3 = 1,
    LN3 = 2,
    EXP = 3
};

// Structure to hold model coefficients and type
struct ModelParams {
    // LogLogistic (3-params) model
    double a = NaN;
    double b = NaN;
    double c = NaN;
    //  LogNormal (3-params) model
    double gamma = NaN;  // Gamma
    double mu = NaN;     // Mean
    double sigma = NaN;  // Standard deviation
    // Exponential model
    double alpha = NaN;  // location of the distribution 
    double beta = NaN;   // 1/beta mean of the distribution
};

struct GoF { // Goodness of fit
    double stat = Inf;  // Statistic
    double thresh = NaN;  // Threshold
};

// Structure to represent a statistical model
struct Model {
    bool defined = false;  // Flag to indicate if the model is defined
    ModelType type = ModelType::None;  // Type of the model (e.g., "LL3", "LN3", "EXP")
    ModelParams params;  // Coefficients of the model
    GoF gof;  // Goodness of fit
};

struct State {
    vector<double> samples = {};  // samples
    Model model;
};

#endif // MODEL_H