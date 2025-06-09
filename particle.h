#pragma once

#include <cassert>              // Include the C++ assert library for runtime assertions.
#include <vector>               // Include the C++ vector library for using std::vector<>.


// Declaration of the Particle class. Objects of this class contain the position vector of the particle
// and its number of dimensions.
//
// Private variables:    unsigned int m_numberOfDimensions - number of dimensions of the particle;
//                       std::vector<double> m_position - position vector of the particle.
//
// Constructor: Particle(const std::vector<double>& position).
//
// Functions that objects of this class can use:    adjustPosition(...);
//                                                  getPosition();
//                                                  getNumberOfDimensions().
class Particle
{
    public:

        //Constructor of the Particle class. It takes a const std::vector<double>& position as its argument,
        //sets the private variable m_position equal to it, and saves the size of this vector to the private variable
        //m_numberOfDimensions.
        Particle(const std::vector<double>& position);


        // Function that adjusts the position of the particle.
        //
        // Input:   double change - the adjustment value;
        //          unsigned int dimension - specifies which coordinate is changed.
        //
        // Output: void - no return value.
        void adjustPosition(double change, unsigned int dimension);


        // Helper function that provides read access to the position vector of the particle.
        //
        // Input:  void.
        //
        // Output: std::vector<double> m_position - the position vector of the particle.
        std::vector<double> &getPosition() { return m_position; }


        // Helper function that provides read access to the number of dimensions of the particle.
        //
        // Input:  void;
        //
        // Output: unsigned int - the number of dimensions of the particle.
        unsigned int getNumberOfDimensions() { return m_numberOfDimensions; }

    private:
        unsigned int m_numberOfDimensions = 0;                    // Number of dimensions of the particle
        std::vector<double> m_position = std::vector<double>();   // Position vector of the particle.
};