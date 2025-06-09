#include <iostream>    // Include the C++ input-output stream library.

#include "initialstate.h"   // Include "initialstate" header file with declarations.



// Function that sets initial particle positions randomly according to a uniform distribution in range [0, 1).
//
// Input:     unsigned int numberOfDimensions - number of dimensions;
//            unsigned int numberOfParticles - number of particles;
//            Random& rng - referance to Random Number Generator.
//
// Output:    std::vector<std::unique_ptr<Particle>> - vector of size numberOfParticles
//                                                     where each element is a unique pointer to a Particle object 
//                                                     containing the coordinates of each particle in each dimension.
std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    // Create a vector that contains unique pointers to Particle class objects.
    auto particles = std::vector<std::unique_ptr<Particle>>();

    // Create a vector that will contain the coordinates of a particle.
    std::vector<double> position = std::vector<double>();

    for (unsigned int i=0; i < numberOfParticles; i++)
    {
        for (unsigned int j=0; j < numberOfDimensions; j++)
        {    
            position.push_back(rng.nextDouble()); // Assign a uniform random number as the j-th coordinate.
        }

        // Assign the position as the coordinate vector of particle i.
        particles.push_back(std::make_unique<Particle>(position));

        position.clear(); // Delete all elements of position vector.
    }

    return particles;
}