#ifndef BASICS_H
#define BASICS_H

#include <cfloat>
#include <cmath>
#include <random>
#include <thread>
#include <vector>
#include <atomic>
#include <mutex>  // Add this include for std::mutex
#include <future> // Add this include for std::future
#include <omp.h>
#include <iostream>
#include <chrono>


std::mutex fileMutex;  // Declare a mutex for file access


constexpr double EPSILON = 1e-6;
constexpr double epsilon = 1e-6;

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Function to generate a random double between min and max
double random_double(double min, double max) {
    static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator; // Mersenne Twister engine
    return distribution(generator);
}

inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}


#endif // BASICS_H
