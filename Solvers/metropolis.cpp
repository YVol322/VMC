#include "metropolis.h"     // Include "metropolis" header file with declarations.


Metropolis::Metropolis(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool Metropolis::step(
            double stepLength,
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        )
{
    int n_particles = particles.size();
    int n_dims = (particles.back()) -> getNumberOfDimensions();

    // Initialize variables to store the wave function values for the old and new configurations.
    double WFold, WFnew;

    // Initialize a 2D vector to store random numbers used to perform particle moves.
    std::vector<std::vector<double>> rand_doubles(n_particles, std::vector<double>(n_dims));

    // Compute the wave function for the initial configuration (old state).
    WFold = waveFunction.evaluate(particles);

    for(int i = 0; i < n_particles; i++)
    {
        for(int j = 0; j < n_dims; j++)
        {
            // Generate a random uniform number between 0 and 1 and save it to the vector.
            rand_doubles[i][j] = m_rng->nextDouble();

            // Move the i-th particle in the j-th dimension by a random amount, based on stepLength and random number.
            (*(particles.at(i))).adjustPosition(stepLength * (rand_doubles[i][j] - 0.5), j);
        }
    }

    // Compute the wave function for the new configuration (after the particle moves).
    WFnew = waveFunction.evaluate(particles);

    // Calculate the acceptance probability based on the Metropolis criterion.
    double accept = std::min(1.0, (WFnew * WFnew) / (WFold * WFold));

    // If a random number is greater than the acceptance probability, reject the step
    // and revert the particle positions back to their original state.
    if(m_rng->nextDouble() > accept)
    {
        // Revert the particle positions using the stored random numbers.
        for(int i = 0; i < n_particles; i++)
        {
            for(int j = 0; j < n_dims; j++)
            {
                (*(particles.at(i))).adjustPosition(-stepLength * (rand_doubles[i][j] - 0.5), j);
            }
        }
        return false; // Step rejected, return false.
    }
    else
    {
        return true; // Step accepted, return true.
    }
}