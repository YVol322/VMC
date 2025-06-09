#pragma once

#include <random>   // Include the C++ random library for generating random numbers.

class Random
{
    private:
        std::mt19937_64 m_engine;  // Mersenne Twister random number generator engine.

    public:
        // Default constructor. Initializes the random number generator with a random seed.
        Random()
        {
            std::random_device rd;      // Obtain a random seed from the system.
            m_engine = std::mt19937_64(rd());  // Initialize the engine with the random seed.
        }



        // Constructor that allows setting a custom seed for the random number generator.
        Random(int seed)
        {
            m_engine = std::mt19937_64(seed);  // Initialize the engine with a user-defined seed.
        }



        // Generates a random integer in the range [lowerLimit, upperLimit].
        // 
        // Input:  const int& lowerLimit - the lower bound of the random number range.
        //         const int& upperLimit - the upper bound of the random number range.
        // 
        // Output: Random integer in the specified range.
        int nextInt(const int &lowerLimit, const int &upperLimit)
        {
            std::uniform_int_distribution<int> dist(lowerLimit, upperLimit);  // Create distribution for the range.
            return dist(m_engine);  // Generate and return the random number.
        }



        // Generates a random integer in the range [0, upperLimit].
        //
        // Input:  const int& upperLimit - the upper bound of the random number range.
        //
        // Output: Random integer in the range [0, upperLimit].
        int nextInt(const int &upperLimit)
        {
            std::uniform_int_distribution<int> dist(0, upperLimit);  // Create distribution from 0 to upperLimit.
            return dist(m_engine);  // Generate and return the random number.
        }



        // Generates a random floating-point number in the range [0, 1).
        //
        // Output: Random floating-point number in the range [0, 1).
        double nextDouble()
        {
            std::uniform_real_distribution<double> dist(0, 1);  // Create distribution for the range [0, 1).
            return dist(m_engine);  // Generate and return the random number.
        }



        // Generates a random number from a Gaussian (normal) distribution.
        // 
        // Input:  const double& mean - the mean (average) of the Gaussian distribution.
        //         const double& standardDeviation - the standard deviation of the Gaussian distribution.
        //
        // Output: Random floating-point number generated from the Gaussian distribution with the specified mean and standard deviation.
        double nextGaussian(const double &mean, const double &standardDeviation)
        {
            std::normal_distribution<double> dist(mean, standardDeviation);  // Create normal distribution.
            return dist(m_engine);  // Generate and return the random number.
        }
};
