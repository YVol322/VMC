#pragma once

#include "montecarlo.h"     // Include "montecarlo" header file with declarations.


// Declaration of the Metropolis class. This class is a subclass of the MonteCarlo class.
// It implements the Metropolis algorithm for Monte Carlo simulations.
//
// Constructor: Metropolis(std::unique_ptr<class Random> rng) - initializes the random number generator for the simulation.
// Function: step(...) - performs a single step of the Metropolis algorithm, generating a new configuration and 
//                       calculating the acceptance probability.
class Metropolis : public MonteCarlo
{
    public:
        // Constructor of the Metropolis class. It initializes the Metropolis simulation with a random number generator.
        //
        // Input:   std::unique_ptr<class Random> rng - a unique pointer to a random number generator instance.
        // Output:  void - no return value.
        Metropolis(std::unique_ptr<class Random> rng);


        // Performs a single step of the Metropolis algorithm.
        //
        // Input:   double stepLength - the length of the step to be taken in the simulation.
        //          class WaveFunction& waveFunction - the wave function object used in the calculation.
        //          std::vector<std::unique_ptr<class Particle>>& particles - a vector of particles in the system.
        //
        // Output:  bool - returns true if the step is accepted, false otherwise.
        bool step(
                double stepLength,
                class WaveFunction& waveFunction,
                std::vector<std::unique_ptr<class Particle>>& particles);
};