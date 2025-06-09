#pragma once

#include <cmath>            // Include the C++ math library for mathematical functions.
#include <Eigen/Dense>      // Include the Eigen library for matrix and vector operations.
#include <iostream>         // Include the C++ input-output stream library.

#include "wavefunction.h"   // Include "wavefunction" header file with declarations.


// Declaration of the Fermion class. This class is a subclass of the WaveFunction class.
// Its purpose is to implement wave function components, gradients, Laplacians, and various Jastrow forms
// for the fermionic system, along with functions to compute quantum forces and evaluate the wave function.
//
// Functions include multiple Psi functions for different components, gradient functions, Laplacians,
// and Jastrow factors, as well as the standard quantum force and local energy computations.
class Fermion : public WaveFunction
{
    public:
        // Constructor of the Fermion class. Initializes the parameters for the fermion system.
        //
        // Input:   double alpha - variational parameter alpha;
        //          std::vector<double> beta - vector of variational parameters for Jastrow or Pade-Jastrow;
        //          int mode - 0 for Jastrow ansatz, 1 for Pade-Jastrow ansatz;
        //          int n_particles - number of particles in the system;
        //          double omega - angular frequency of the system.
        Fermion(double alpha, std::vector<double> beta, int mode, int n_particles, double omega);
    

        // Computes the Jastrow factor for the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed Jastrow factor for the given particle configuration.
        double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);


        // Computes the Pade-Jastrow factor for the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed Pade-Jastrow factor for the given particle configuration.
        double PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles);


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


        // Computes the Slater determinant or its gradient component or Laplacian depending on the derivative order.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int particles_set - set of particles: 0 for spin-up, 1 for spin-down;
        //          int row_changed - row index of the particle that was changed;
        //          int der_order - order of the derivatives to compute: 0 for Psi, 1 for gradient, 2 for Laplacian;
        //          int grad_comp - component of the gradient to compute (if applicable).
        //
        // Output:  double - the computed value based on the requested derivative order (Slater determinant, gradient, or Laplacian).
        double SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp);


        // Computes the gradient of the Jastrow factor for a given particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of the Jastrow factor.
        //
        // Output:  std::vector<double> - the gradient of the Jastrow factor for the given particle.
        std::vector<double> GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the gradient of the Pade-Jastrow factor for a given particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_inx - index of the particle for which to compute the gradient of the Pade-Jastrow factor.
        //
        // Output:  std::vector<double> - the gradient of the Pade-Jastrow factor for the given particle.
        std::vector<double> GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);


        // Computes the Laplacian of the Jastrow factor for a given particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the Laplacian of the Jastrow factor.
        //
        // Output:  double - the Laplacian of the Jastrow factor for the given particle.
        double LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // Computes the Laplacian of the Pade-Jastrow factor for a given particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_inx - index of the particle for which to compute the Laplacian of the Pade-Jastrow factor.
        //
        // Output:  double - the Laplacian of the Pade-Jastrow factor for the given particle.
        double LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);



        // Computes the sum of Laplacians of Psi_T divided by Psi_T for all particles in the system using the Fermion ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
        double LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles);


        // Computes the gradient of Psi_T multiplied by the gradient of the Jastrow factor (or the Pade-Jastrow factor) 
        // over Psi_T for a given particle.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          double part_idx - index of the particle for which to compute the gradient of Psi_T multiplied by
        //                            gradient of J (or PJ) over Psi_T.
        //
        // Output:  double - the computed value of the gradient of Psi_T multiplied by gradient of J (or PJ) over Psi_T.
        double GradiPsiTGradiJOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);


        // WaveFunction class function that evaluates the wave function for the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed value of the wave function for the given particle configuration.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
        // for the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the quantum force for a specific particle in the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int i - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the quantum force for the specified particle.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int i);
};