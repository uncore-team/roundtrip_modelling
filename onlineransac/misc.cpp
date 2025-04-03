#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

/**
 * Writes numerical data from a vector into a file.
 *
 * @param filename Path to the output file
 * @param data Vector of double values to be written to the file
 *
 * Implementation details:
 * - Opens file in text mode
 * - Writes space-separated double values
 * - Validates file opening
 *
 * @throws runtime_error if file cannot be opened
 */
void write_data(const string& filename, const vector<double>& data) {
    ofstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file");
    }
    for (const double& value : data) {
        file << value << endl;
    }
    file.close();
}

/**
 * Reads numerical data from a file into a vector.
 *
 * @param filename Path to the input file
 * @return Vector of double values read from file
 *
 * Implementation details:
 * - Opens file in text mode
 * - Reads space-separated double values
 * - Validates file opening
 *
 * @throws runtime_error if file cannot be opened
 */
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
