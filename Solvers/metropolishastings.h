#pragma once

#include "montecarlo.h"     // Include "montecarlo" header file with declarations.


// Declaration of the MetropolisHastings class, which inherits from the MonteCarlo base class.
// This class implements the Metropolis-Hastings algorithm for Monte Carlo simulations, which includes
// a method to perform a single step of the simulation by proposing a particle move and accepting or rejecting it
// based on the acceptance probability calculated using the wave function.
class MetropolisHastings : public MonteCarlo
{
    public:
        // Constructor for the MetropolisHastings class. It initializes the random number generator (rng).
        // Input:   std::unique_ptr<class Random> rng - a random number generator to be used in the simulation.
        MetropolisHastings(std::unique_ptr<class Random> rng);

        // Method to perform a single step in the Metropolis-Hastings simulation.
        // The algorithm proposes a move for each particle and accepts or rejects the move based on the 
        // Metropolis-Hastings criterion.
        //
        // Input:   double stepLength - the magnitude of the proposed step for particle movement.
        //          class WaveFunction& waveFunction - reference to the wave function used to calculate the 
        //                                             acceptance probability.
        //          std::vector<std::unique_ptr<class Particle>>& particles - list of particles to be moved.
        //
        // Output:  bool - true if the step is accepted, false if the step is rejected.
        bool step(
                double stepLength,
                class WaveFunction& waveFunction,
                std::vector<std::unique_ptr<class Particle>>& particles);
};