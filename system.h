#pragma once

#include <cassert>              // Include the C++ assert library for runtime assertions.
#include <memory>               // Include the C++ memory library for using std::unique_ptr.
#include <vector>               // Include the C++ vector library for using std::vector<>.

#include "sampler.h"                            // Include "sampler" header file with declarations.
#include "particle.h"                           // Include "particle" header file with declarations.
#include "WaveFunctions/wavefunction.h"         // Include "wavefunction" header file with declarations.
#include "Hamiltonians/hamiltonian.h"           // Include "hamiltonian" header file with declarations.
#include "InitialStates/initialstate.h"         // Include "initialstate" header file with declarations.
#include "Solvers/montecarlo.h"                 // Include "montecarlo" header file with declarations.


// Declaration of the System class. This class represents a system of particles that includes 
// the Hamiltonian, wave function, Monte Carlo solver, and particle configurations.
//
// The class provides functionality to perform equilibration steps, run Metropolis steps, 
// compute local energy, and other particle-system specific operations.
//
// Private variables:    unsigned int m_numberOfParticles - number of particles in the system;
//                       unsigned int m_numberOfDimensions - number of dimensions of the system;
//                       std::unique_ptr<class Hamiltonian> m_hamiltonian - Hamiltonian of the system;
//                       std::unique_ptr<class WaveFunction> m_waveFunction - wave function for the system;
//                       std::unique_ptr<class MonteCarlo> m_solver - solver for the Monte Carlo simulation;
//                       std::vector<std::unique_ptr<class Particle>> m_particles - list of particles in the system.
//
// Constructor: System(...).
//
// Functions:   runEquilibrationSteps(...), 
//              runMetropolisSteps(...), 
//              computeLocalEnergy(), 
//              computerij(), 
//              computer2(), 
//              getWaveFunctionParameters(), 
//              etc.
class System
{
        public:
        // Constructor for the System class. Initializes the system with the given Hamiltonian, wave function,
        // Monte Carlo solver, and list of particles.
        //
        // Input:   std::unique_ptr<class Hamiltonian> hamiltonian - Hamiltonian of the system;
        //          std::unique_ptr<class WaveFunction> waveFunction - wave function of the system;
        //          std::unique_ptr<class MonteCarlo> solver - solver for the Monte Carlo simulation;
        //          std::vector<std::unique_ptr<class Particle>> particles - list of particles.
        System(
                std::unique_ptr<class Hamiltonian> hamiltonian,
                std::unique_ptr<class WaveFunction> waveFunction,
                std::unique_ptr<class MonteCarlo> solver,
                std::vector<std::unique_ptr<class Particle>> particles);


        // Performs equilibration steps for the system, adjusting the state based on the given step length
        // and number of equilibration steps.
        //
        // Input:   double stepLength - step size for each Monte Carlo step;
        //          unsigned int numberOfEquilibrationSteps - number of equilibration steps to perform.
        //
        // Output:  unsigned int - number of equilibration steps performed.
        unsigned int runEquilibrationSteps(
                double stepLength,
                unsigned int numberOfEquilibrationSteps);


        // Runs Metropolis steps for the system and returns a sampler object for collecting data during the simulation.
        //
        // Input:   double stepLength - step size for each Monte Carlo step;
        //          unsigned int numberOfMetropolisSteps - number of Metropolis steps to perform.
        //
        // Output:  std::unique_ptr<class Sampler> - a sampler for collecting data during the simulation.
        std::unique_ptr<class Sampler> runMetropolisSteps(
                double stepLength,
                unsigned int numberOfMetropolisSteps);


        // Computes the local energy of the system, which is required for variational Monte Carlo simulations.
        //
        // Output:  double - computed local energy of the system.
        double computeLocalEnergy();


        // Computes the distance between particles i and j in the system.
        //
        // Input:   int i - index of the first particle;
        //          int j - index of the second particle.
        //
        // Output:  double - distance between particles i and j.
        double computerij(int i, int j);


        // Computes the sum of the square of the distances (r^2) for all particles in the system.
        //
        // Output:  double - sum of r^2 for all particles in the system.
        double computer2();


        // Helper function that provides access to the current wave function parameters.
        //
        // Output:  const std::vector<double>& - the current wave function parameters.
        const std::vector<double>& getWaveFunctionParameters();

        private:
            unsigned int m_numberOfParticles = 0;
            unsigned int m_numberOfDimensions = 0;

            std::unique_ptr<class Hamiltonian> m_hamiltonian;    // Hamiltonian for the system.
            std::unique_ptr<class WaveFunction> m_waveFunction;  // Wave function for the system.
            std::unique_ptr<class MonteCarlo> m_solver;         // Monte Carlo solver for the system.
            std::vector<std::unique_ptr<class Particle>> m_particles;  // List of particles in the system.
};