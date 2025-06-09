#pragma once
#include <cmath>                // Include the C++ math library for mathematical functions.
#include <memory>               // Include the C++ memory library for using std::unique_ptr.
#include <vector>               // Include the C++ vector library for using std::vector<>.

#include "system.h"                         // Include "system" header file with declarations.
#include "sampler.h"                        // Include "sampler" header file with declarations.
#include "particle.h"                       // Include "particle" header file with declarations.
#include "Hamiltonians/hamiltonian.h"       // Include "hamiltonian" header file with declarations.
#include "WaveFunctions/wavefunction.h"     // Include "wavefunction" header file with declarations.


// Declaration of the Sampler class. This class is responsible for handling the sampling process in a 
// Variational Monte Carlo (VMC) simulation, including the computation of averages for observables 
// like local energy, Jastrow factor, and others. It also manages the number of Metropolis steps and 
// stores the results for optimization of variational parameters like alpha, beta, and beta_ij.
//
//
// Constructor:        Sampler(unsigned int numberOfParticles, unsigned int numberOfDimensions, 
//                       double stepLength, unsigned int numberOfMetropolisSteps);
// 
// Functions:          sample() - performs the sampling of observables;
//                     computeAverages() - computes the averages of sampled variables;
//                     printOutputToTerminal() - prints the results to the terminal;
//                     setEnergy() - sets the energy for the VMC simulation;
//                     setTime() - sets the simulation time;
//                     getEnergy(), getO1alpha(), getO2alpha(), getO1Jastow(), getO2Jastow() - 
//                     accessors for the averages of various observables.
class Sampler
{
    public:
        // Constructor to initialize the sampler with necessary parameters.
        //
        // Input:   unsigned int numberOfParticles - the number of particles in the system;
        //          unsigned int numberOfDimensions - the number of dimensions in the system;
        //          double stepLength - the step length used in the Metropolis algorithm;
        //          unsigned int numberOfMetropolisSteps - total number of Metropolis steps to sample.
        Sampler(
            unsigned int numberOfParticles,
            unsigned int numberOfDimensions,
            double stepLength,
            unsigned int numberOfMetropolisSteps);


        // Samples all observables, including the local energy, O1, and O2 for alpha, beta, and beta_ij optimization.
        //
        // Input:   bool acceptedStep - indicates if the step was accepted in the Metropolis algorithm;
        //          class System* system - pointer to the system object that contains the simulation state.
        void sample(bool acceptedStep, class System* system);


        // Prints the output to the terminal, typically the results of the VMC simulation.
        //
        // Input:   class System& system - reference to the system object to extract the necessary information.
        void printOutputToTerminal(class System& system);


        // Computes the averages of all sampled variables after the sampling is complete.
        void computeAverages();


        // Sets the VMC energy to be used in output (necessary for parallel algorithms).
        //
        // Input:   double en - the value of the energy to set.
        void setEnergy(double en);
        

        // Sets the VMC time to be used in output (necessary for parallel algorithms).
        //
        // Input:   double t - the runtime of the VMC simulation.
        void setTime(double t);

        // Helper functions to collect mean values of various observables:
        // Energy, O1, and O2 for various optimizations.
        double getEnergy() { return m_energy; }        // Returns the mean local energy.
        double getO1alpha() { return m_O1alpha; }      // Returns the mean value of O1 for alpha optimization.
        double getO2alpha() { return m_O2alpha; }      // Returns the mean value of O2 for alpha optimization.
        std::vector<double> getO1Jastow() { return m_O1Jastrow; }  // Returns the mean O1 for Jastrow optimization.
        std::vector<double> getO2Jastow() { return m_O2Jastrow; }  // Returns the mean O2 for Jastrow optimization.
        double getO1Pade() { return m_O1Pade; }        // Returns the mean O1 for Pade-Jastrow optimization.
        double getO2Pade() { return m_O2Pade; }        // Returns the mean O2 for Pade-Jastrow optimization.

    private:
        unsigned int m_stepNumber = 0;                 // Current iteration number in the VMC simulation.
        unsigned int m_numberOfMetropolisSteps = 0;    // Total number of Metropolis steps.
        unsigned int m_numberOfParticles = 0;          // Number of particles in the system.
        unsigned int m_numberOfDimensions = 0;         // Number of dimensions of the system.
        unsigned int m_numberOfAcceptedSteps = 0;      // Number of accepted Metropolis steps.

        double m_stepLength = 0;                       // Step length used in the Metropolis algorithm.
        double m_time = 0;                             // VMC algorithm runtime.
        int m_nPairs = 0;                              // Number of Jastrow factor variational parameters (beta_ij).

        double m_energy = 0;                           // Stores the mean local energy for the system.
        double m_cumulativeEnergy = 0;                 // Cumulative sum of the local energies for averaging.

        double m_O1alpha = 0;                          // Stores the mean O1 for alpha optimization.
        double m_cumulativeO1alpha = 0;                // Cumulative sum for O1 alpha optimization.

        double m_O2alpha = 0;                          // Stores the mean O2 for alpha optimization.
        double m_cumulativeO2alpha = 0;                // Cumulative sum for O2 alpha optimization.

        std::vector<double> m_O1Jastrow;               // Stores the mean O1 for Jastrow factor optimization.
        std::vector<double> m_cumulativeO1Jastrow;     // Cumulative sum for O1 Jastrow optimization.

        std::vector<double> m_O2Jastrow;               // Stores the mean O2 for Jastrow factor optimization.
        std::vector<double> m_cumulativeO2Jastrow;     // Cumulative sum for O2 Jastrow optimization.

        double m_O1Pade;                               // Stores the mean O1 for Pade-Jastrow optimization.
        double m_cumulativeO1Pade = 0;                 // Cumulative sum for O1 Pade-Jastrow optimization.

        double m_O2Pade;                               // Stores the mean O2 for Pade-Jastrow optimization.
        double m_cumulativeO2Pade = 0;                 // Cumulative sum for O2 Pade-Jastrow optimization.
};