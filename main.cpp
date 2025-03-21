#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "OnlineRANSAC.h"
#include "ExponentialEstimator.h"
// #if _PYTHON
//     #include "matplotlibcpp.h"
//     namespace plt = matplotlibcpp;
// #endif

using namespace std;

/**
 * Main program entry point.
 * Tests Online RANSAC algorithm with different configurations.
 * 
 * Implementation details:
 * - Reads RTT data from external file
 * - Tests combinations of:
 *   * model_preserving (0-2)
 *   * sample_sliding (0-1)
 *   * data_preserving (0-1)
 * - Uses LN3 model type
 * - Prints results for each configuration
 * 
 * @return 0 on successful execution
 */
int main() {

    // Read RTTs from a file
    const string filename = "matlab/rtts.txt";
    const vector<double> samples = read_data(filename);

    // List the estimators
    const vector<ModelType> model_types = {ModelType::LL3};

    const int len = samples.size();
    for(int model_preserving = 0; model_preserving < 3; ++model_preserving) {
        for (int sample_sliding = 0; sample_sliding < 2; ++sample_sliding) {
            for (int data_preserving = 0; data_preserving < 2; ++data_preserving) {

                // Initilize the Online RANSAC algorithm
                OnlineRANSAC onlineRANSAC(
                    20, (unsigned)Inf,
                    model_preserving, sample_sliding, data_preserving,
                    model_types
                );

                for (int f = 0; f < len; ++f) {

                    auto t1 = chrono::high_resolution_clock::now();
                    auto exitbranch = onlineRANSAC.update(samples[f]);
                    auto t2 = chrono::high_resolution_clock::now();
                    chrono::duration<double> ct = t2 - t1;
 //                   cout << "#" << (f+1) << ", exitbranch[" << exitbranch << "], time[" << ct.count() << "]" << endl;
                }

                // print the final model if it exists
                cout << endl << "model_preserving[" << model_preserving << "] "
                    << "sample_sliding[" << sample_sliding << "] "
                    << "data_preserving[" << data_preserving << "]";
                onlineRANSAC.print_model();
            }
        }
    }

// #ifdef _PYTHON
//     // Initialize Python environment
//     if (!plt::PythonEnvironment::initialize()) {
//         std::cerr << "Error initializing Python" << std::endl;
//         return 1;
//     }
    
//     // Prepare data.
//     int n = 5000;
//     std::vector<double> x(n), y(n), z(n), w(n,2);
//     for(int i=0; i<n; ++i) {
//         x.at(i) = i*i;
//         y.at(i) = sin(2*M_PI*i/360.0);
//         z.at(i) = log(i);
//     }

//     // Set the size of output image to 1200x780 pixels
//     plt::figure_size(1200, 780);
//     // Plot line from given x and y data. Color is selected automatically.
//     plt::plot(x, y);
//     // Plot a red dashed line from given x and y data.
//     plt::plot(x, w,"r--");
//     // Plot a line whose name will show up as "log(x)" in the legend.
//     plt::named_plot("[FAKE] log(x)", x, z);
//     // Set x-axis to interval [0,1000000]
//     plt::xlim(0, 1000*1000);
//     // Add graph title
//     plt::title("[FAKE] Sample figure");
//     // Enable legend.
//     plt::legend();
//     // Show/Save the image (file format is determined by the extension)
//     plt::show(); //    plt::save("./basic.png");

//     plt::PythonEnvironment::finalize();
// #endif

    return 0;
}