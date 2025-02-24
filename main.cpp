#include <fstream>
#include <stdexcept>
#include <chrono>
#include <iostream>
#include <vector>
#include <string>

#include "OnlineRANSAC.h"
#include "ExponentialEstimator.h"

using namespace std;

vector<double> read_data(const string& filename) {
    vector<double> samples;
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file");
    }
    double value;
    while (file >> value) {
        samples.push_back(value);
    }
    file.close();
    return samples;
}

string model_to_string(ModelType model_type) {
    switch(model_type) {
        case ModelType::LL3: return "LL3"; break;
        case ModelType::LN3: return "LN3"; break;
        case ModelType::EXP: return "EXP"; break;
        case ModelType::None: return "None"; break;
    }
    return "Unknown";
}

void print_model(const Model& model) {
    if (!model.defined) {
        cout << "Model is not defined." << endl;
        return;
    }

    cout << "Model Type: " << model_to_string(model.type) << endl;
    cout << "Coefficients:" << endl;
    cout << "a: " << model.params.a << endl;
    cout << "b: " << model.params.b << endl;
    cout << "c: " << model.params.c << endl;
    cout << "alpha: " << model.params.alpha << endl;
    cout << "beta: " << model.params.beta << endl;
    cout << "gamma: " << model.params.gamma << endl;
    cout << "mu: " << model.params.mu << endl;
    cout << "sigma: " << model.params.sigma << endl;
}

int main() {

    // Read RTTs from a file
    string filename = "rtts.txt";
    vector<double> samples = read_data(filename);

    // Prepare the estimators
    vector<ModelType> model_types = {ModelType::EXP};

    // Initilize the Online RANSAC algorithm
    OnlineRANSAC onlineRANSAC(20, (unsigned)INFINITY, false, false, false, model_types);

    int len = samples.size();

    for (int f = 0; f < len; ++f) {

//        auto t1 = chrono::high_resolution_clock::now();
        auto exitbranch = onlineRANSAC.update(samples[f]);
        cout << "#" << (f+1) << "exitbranch[" << exitbranch << "]" << endl;
//        auto t2 = chrono::high_resolution_clock::now();
//        chrono::duration<double> ct = t2 - t1;
    }

    // Get the final model
    Model model = onlineRANSAC.get_model();

    print_model(model);

    return 0;
}