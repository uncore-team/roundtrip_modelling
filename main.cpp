#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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
        cout << "\nModel Type: NOT defined." << endl;
        return;
    }

    cout << "\nModel Type: " << model_to_string(model.type);
    cout << fixed << setprecision(6);  // Configurar formato para números flotantes
    
    switch(model.type) {
        case ModelType::LL3:
            cout << "\n\ta: " << model.params.a
                 << "\n\tb: " << model.params.b
                 << "\n\tc: " << model.params.c;
            break;
        case ModelType::LN3:
            cout << "\n\tgamma: " << model.params.gamma
                 << "\n\tmu: " << model.params.mu
                 << "\n\tsigma: " << model.params.sigma;
            break;
        case ModelType::EXP:
            cout << "\n\talpha: " << model.params.alpha
                 << "\n\tbeta: " << model.params.beta;
            break;
        case ModelType::None:
            cout << "Model is not defined.";
            break;
    }
    cout << endl;
}

int main() {

    // Read RTTs from a file
    string filename = "matlab/rtts.txt";
    vector<double> samples = read_data(filename);

    // List the estimators
    vector<ModelType> model_types = {ModelType::LL3};

    int len = samples.size();
    for(int model_preserving = 0; model_preserving < 3; ++model_preserving) {
        for (int sample_sliding = 0; sample_sliding < 2; ++sample_sliding) {
            for (int data_preserving = 0; data_preserving < 2; ++data_preserving) {

                cout << endl << "model_preserving[" << model_preserving << "] "
                    << "sample_sliding[" << sample_sliding << "] "
                    << "data_preserving[" << data_preserving << "]";

                // Initilize the Online RANSAC algorithm
                OnlineRANSAC onlineRANSAC(
                    20, (unsigned)Inf,
                    model_preserving, sample_sliding, data_preserving,
                    model_types
                );
                for (int f = 0; f < len; ++f) {

                    // auto t1 = chrono::high_resolution_clock::now();
                    auto exitbranch = onlineRANSAC.update(samples[f]);
                    // auto t2 = chrono::high_resolution_clock::now();
                    // chrono::duration<double> ct = t2 - t1;
                    // cout << "#" << (f+1) << ", exitbranch[" << exitbranch << "], time[" << ct.count() << "]" << endl;
                }

                // Get the final model
                print_model(onlineRANSAC.get_model());

//                return 0; // just for testing the first iteration
            }
        }
    }

    return 0;
}