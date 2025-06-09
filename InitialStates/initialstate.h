#pragma once

#include <cassert>              // Include the C++ assert library for runtime assertions.
#include <memory>               // Include the C++ memory library for using std::unique_ptr.
#include <vector>               // Include the C++ vector library for using std::vector<>.

#include "../particle.h"        // Include "particle" header file with declarations.
#include "../Math/random.h"     // Include "random" header file with declarations.



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
            Random& randomEngine
        );