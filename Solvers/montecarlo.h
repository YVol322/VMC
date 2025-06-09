#pragma once

#include <memory>               // Include the C++ memory library for std::unique_ptr.
#include <vector>               // Include the C++ vector library for using std::vector<>.

#include "../Math/random.h"                     // Include "random" header file with declarations.
#include "../WaveFunctions/wavefunction.h"      // Include "wavefunction" header file with declarations.
#include "../particle.h"                        // Include "particle" header file with declarations.


// Declaration of the MonteCarlo class. This class provides a base class for implementing different Monte Carlo simulation algorithms.
//
// Private variables:    std::unique_ptr<class Random> m_rng - a unique pointer to a random number generator instance.
//
// Constructor: MonteCarlo(std::unique_ptr<class Random> rng) - initializes the Monte Carlo simulation with a random number generator.
//
// Virtual functions that derived classes need to implement:    step(double stepLength, class WaveFunction& waveFunction, 
//                                                              std::vector<std::unique_ptr<class Particle>>& particles) - 
//                                                              performs a single step of the simulation (pure virtual function).
class MonteCarlo
{
    public:

        // Constructor of the MonteCarlo class. It initializes the Monte Carlo simulation with a random number generator.
        //
        // Input:   std::unique_ptr<class Random> rng - a unique pointer to a random number generator instance.
        //
        // Output:  void - no return value.
        MonteCarlo(std::unique_ptr<class Random> rng);


        // Virtual destructor for the MonteCarlo class. The virtual keyword ensures that the correct derived class destructor is called.
        // The destructor is marked default, meaning the compiler will automatically generate it.
        virtual ~MonteCarlo() = default;


        // Pure virtual function 'step' that performs a single Monte Carlo step. This function must be implemented by derived classes.
        //
        // Input:   double stepLength - the step size or distance to move the particles.
        //          class WaveFunction& waveFunction - a reference to the wavefunction used to compute energy or accept/reject the step.
        //          std::vector<std::unique_ptr<class Particle>>& particles - a list of particles in the system.
        //
        // Output:  bool - returns true if the step was accepted, false otherwise.
        virtual bool step(
                double stepLength,
                class WaveFunction& waveFunction,
                std::vector<std::unique_ptr<class Particle>>& particles) = 0;

    protected:
        // Member variable to hold the random number generator instance.
        // The random number generator is used by subclasses to generate random numbers for the simulation.
        std::unique_ptr<class Random> m_rng;
};