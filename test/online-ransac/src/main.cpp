#include <iostream>
#include "online_ransac.h"

int main() {
    // Initialize the OnlineRANSAC class
    OnlineRANSAC ransac;

    // Example data for estimation
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};

    // Estimate the exponential distribution
    ransac.estimate(data);

    // Output the estimated parameters
    ransac.printParameters();

    return 0;
}