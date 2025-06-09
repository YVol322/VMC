#pragma once

#include <cmath>            // Include the C++ math library for mathematical functions.
#include <iostream>         // Include the C++ input-output stream library.

#include "wavefunction.h"   // Include "wavefunction" header file with declarations.


// Declaration of the Boson class. This class is a subclass of the WaveFunction class.
// Its purpose is to implement the evaluate, computeDoubleDerivative, and quantumForce functions for the Bosonic ansatz.
// All other functions are designed to simplify these three core functions.
class Boson : public WaveFunction
{
    public:
        // Constructor of the Boson class. Initializes the parameters for the boson system.
        //
        // Input:   double alpha - variational parameter alpha;
        //          int n_particles - number of particles in the system;
        //          double omega - angular frequency of the system.
        Boson(double alpha, int n_particles, double omega);


        // Boson subclass specific function that computes the trial wave function (Psi_T).
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the value of the trial wave function (Psi_T) for the given particle configuration.
        double PsiT(std::vector<std::unique_ptr<class Particle>>& particles);


        // Boson subclass specific function that computes the sum of Laplacians of Psi_T divided by Psi_T for all particles.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
        double LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that evaluates the wave function for the given particle configuration
        // using the bosonic ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the value of the wave function for the given particle configuration using the bosonic ansatz.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
        // for the bosonic ansatz, which is required for calculating the local energy.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the quantum force for a specific particle in the system
        // using the bosonic ansatz. The quantum force is used in the Metropolis-Hastings sampling algorithm.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int n - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the quantum force for the specified particle.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
};