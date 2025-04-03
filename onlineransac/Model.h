#ifndef MODEL_H
#define MODEL_H

#include <limits>

using namespace std;

#define MIN(a,b) ((a) < (b) ? a : b)
#define MAX(a,b) ((a) > (b) ? a : b)

static const double eps = numeric_limits<double>::epsilon();
static const double Inf = numeric_limits<double>::max();
static const double NaN = numeric_limits<double>::quiet_NaN();
static const size_t OMP_THRESH = 1000;

// Model types
enum ModelType {
    None = 0,
    LL3 = 1,
    LN3 = 2,
    EXP = 3
};

// Structure to hold model coefficients
struct ModelParams {
    // LogLogistic (3-params) model
    double a = NaN;      // location parameter (minimum possible value)
    double b = NaN;      // scale parameter (alpha)
    double c = NaN;      // shape parameter (sigma)
    //  LogNormal (3-params) model
    double gamma = NaN;  // Gamma
    double mu = NaN;     // Mean
    double sigma = NaN;  // Standard deviation
    // Exponential model
    double alpha = NaN;  // location of the distribution
    double beta = NaN;   // 1/beta mean of the distribution
};
//typedef vector<double> ModelParams; TODO (a vector instead of the previous struct)

// Goodness of fit
struct GoF {
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

extern void write_data(const string& filename, const vector<double>& data);
extern vector<double> read_data(const string& filename);


#endif // MODEL_H
