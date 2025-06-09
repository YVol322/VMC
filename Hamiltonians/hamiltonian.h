#pragma once

#include <cassert>                              // Include the C++ assert library for runtime assertions
#include <memory>                               // Include the C++ memory library for using std::unique_ptr
#include <vector>                               // Include the C++ vector library for using std::vector<>

#include "../particle.h"                        // Include "particle" header file with declarations.
#include "../WaveFunctions/wavefunction.h"      // Include "wavefunction" header file with declarations.



// Declaration of the Hamiltonian class. This class is designed to compute the local energy,
// which is sampled during VMC iterations. Hamiltonian parts are constructed as subclasses of the Hamiltonian class.
//
//
// Functions that objects of this class can use:    computeLocalEnergy(...);
//
//
// Virtual destructor: ~Hamiltonian() - virtual destructor for proper cleanup of derived class objects.
class Hamiltonian
{
    public:
        // Virtual destructor of the Hamiltonian class. Ensures proper cleanup of derived class objects.
        virtual ~Hamiltonian() = default;


        // Pure virtual function that computes the local energy for the particle system.
        // This function must be implemented by any derived class.
        //
        // Input:   class WaveFunction& waveFunction - the wavefunction of the system;
        //          std::vector<std::unique_ptr<class Particle>>& particles - the list of particles in the system.
        //
        // Output:  double - the computed local energy.
        virtual double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        ) = 0;
};
