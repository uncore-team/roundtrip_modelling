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
        printf("Model is not defined.\n");
        return;
    }

    printf("\nModel Type: %s", model_to_string(model.type).c_str());
    switch(model.type) {
        case ModelType::LL3:
            printf("\n\ta: %.6f", model.params.a);
            printf("\n\tb: %.6f", model.params.b);
            printf("\n\tc: %.6f", model.params.c);
            break;
        case ModelType::LN3:
            printf("\n\tgamma: %.6f", model.params.gamma);
            printf("\n\tmu: %.6", model.params.mu);
            printf("\n\tsigma: %.6f", model.params.sigma);
            break;
        case ModelType::EXP:
            printf("\n\talpha: %.6f", model.params.alpha);
            printf("\n\tbeta: %.6f", model.params.beta);
            break;
        case ModelType::None:
            printf("Model is not defined.\n");
            break;
    }
}

int main() {

    // Read RTTs from a file
    string filename = "rtts.txt";
    vector<double> samples = read_data(filename);

    // Prepare the estimators
    vector<ModelType> model_types = {ModelType::EXP};

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

            //        auto t1 = chrono::high_resolution_clock::now();
                    auto exitbranch = onlineRANSAC.update(samples[f]);
                    // cout << "#" << (f+1) << " exitbranch[" << exitbranch << "]" << endl;
            //        auto t2 = chrono::high_resolution_clock::now();
            //        chrono::duration<double> ct = t2 - t1;
                }

                // Get the final model
                print_model(onlineRANSAC.get_model());       
            }
        }
    }

    return 0;
}