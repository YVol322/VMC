#pragma once

#include <cmath>            // Include the C++ math library for mathematical functions.
#include <Eigen/Dense>      // Include the Eigen library for matrix and vector operations.
#include <iostream>         // Include the C++ input-output stream library.

#include "wavefunction.h"   // Include "wavefunction" header file with declarations.


// Declaration of the FermionNI class. This class is a subclass of the WaveFunction class.
// Its purpose is to implement wave function components, gradients, and Laplacians for non-interacting
// fermionic system, as well as functions for computing quantum forces and evaluating the wave function.
//
// Functions include multiple Psi functions for different components, gradient functions, Laplacians, 
// and the standard quantum force and local energy computations.
class FermionNI : public WaveFunction
{
    public:
        // Constructor of the FermionNI class. Initializes the parameters for the fermion system.
        //
        // Input:   double alpha - variational parameter alpha;
        //          int n_particles - number of particles in the system;
        //          double omega - angular frequency of the system.
        FermionNI(double alpha, int n_particles, double omega);
    

        // Computes the value of Psi1 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute Psi1.
        //
        // Output:  double - the value of Psi1 for the given particle.
        double Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the value of Psi2 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute Psi2.
        //
        // Output:  double - the value of Psi2 for the given particle.
        double Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the value of Psi3 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute Psi3.
        //
        // Output:  double - the value of Psi3 for the given particle.
        double Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the value of Psi4 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute Psi4.
        //
        // Output:  double - the value of Psi4 for the given particle.
        double Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the value of Psi5 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute Psi5.
        //
        // Output:  double - the value of Psi5 for the given particle.
        double Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the value of Psi6 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute Psi6.
        //
        // Output:  double - the value of Psi6 for the given particle.
        double Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of Psi1 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi1.
        //
        // Output:  std::vector<double> - the gradient of Psi1 for the given particle.
        std::vector<double> GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of Psi2 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi2.
        //
        // Output:  std::vector<double> - the gradient of Psi2 for the given particle.
        std::vector<double> GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of Psi3 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi3.
        //
        // Output:  std::vector<double> - the gradient of Psi3 for the given particle.
        std::vector<double> GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of Psi4 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi4.
        //
        // Output:  std::vector<double> - the gradient of Psi4 for the given particle.
        std::vector<double> GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of Psi5 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi5.
        //
        // Output:  std::vector<double> - the gradient of Psi5 for the given particle.
        std::vector<double> GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of Psi6 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi6.
        //
        // Output:  std::vector<double> - the gradient of Psi6 for the given particle.
        std::vector<double> GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of Psi1 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of Psi1.
        //
        // Output:  double - the Laplacian of Psi1 for the given particle.
        double LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of Psi2 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of Psi2.
        //
        // Output:  double - the Laplacian of Psi2 for the given particle.
        double LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of Psi3 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of Psi3.
        //
        // Output:  double - the Laplacian of Psi3 for the given particle.
        double LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of Psi4 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of Psi4.
        //
        // Output:  double - the Laplacian of Psi4 for the given particle.
        double LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of Psi5 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of Psi5.
        //
        // Output:  double - the Laplacian of Psi5 for the given particle.
        double LapliPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of Psi6 for the specified particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of Psi6.
        //
        // Output:  double - the Laplacian of Psi6 for the given particle.
        double LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the spin-up Slater determinant, its gradient component, or Laplacian depending on the arguments.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int particles_set - number of sets of particles: 0 - spin-up, 1 - spin-down;
        //          int row_changed - row index of the particle that was changed;
        //          int der_order - order of the derivatives to compute;
        //          int grad_comp - the gradient component to compute (if applicable).
        //
        // Output:  double - the computed value based on the requested derivative order (Slater determinant, gradient, or Laplacian).
        double SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp);


        // Computes the sum of Laplacians of Psi_T divided by Psi_T for all particles in the system using the Fermion ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
        double LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the quantum force for a specific particle in the system
        // using the fermionic ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int n - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the quantum force for the specified particle.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);


        // WaveFunction class function that evaluates the wave function for the given particle configuration 
        // using the fermionic ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the value of the wave function for the given particle configuration using the fermionic ansatz.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
        // for the fermionic ansatz, which is required for calculating the local energy.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
};