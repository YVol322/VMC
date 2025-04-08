#include <iostream>

#include "initialstate.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();

    std::vector<double> position = std::vector<double>();

    for (unsigned int i=0; i < numberOfParticles; i++)
    {
        for (unsigned int j=0; j < numberOfDimensions; j++)
        {    
            position.push_back(rng.nextDouble());

            //std::cout << "particle number " << i + 1 << ", dimnsion " << j + 1  << "coordinate " << position.at(j)<< std::endl;
        }

        particles.push_back(std::make_unique<Particle>(position));

        position.clear();
    }

    return particles;
}