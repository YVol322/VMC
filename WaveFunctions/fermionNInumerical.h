#pragma once

#include <cmath>            // Include the C++ math library for mathematical functions.
#include <Eigen/Dense>      // Include the Eigen library for matrix and vector operations.
#include <iostream>         // Include the C++ input-output stream library.

#include "wavefunction.h"   // Include "wavefunction" header file with declarations.

#include <autodiff/reverse/var.hpp>         // Include autodiff reverse-mode library for automatic differentiation.
#include <autodiff/reverse/var/eigen.hpp>   // Include Eigen support for autodiff reverse-mode variables.
using namespace autodiff;
using autodiff::var;
using autodiff::VectorXvar;
using autodiff::MatrixXvar;


// Declaration of the FermionNInumerical class. This class is a subclass of the WaveFunction class.
// It implements the wave function components, gradients, and Laplacians for a non-interacting fermionic system 
// using automatic differentiation. It also provides functions for computing quantum forces and evaluating the wave function.
//
// Functions include multiple psi functions for different components (psi1i, psi2i, ...), gradient functions, 
// Laplacians, and the standard quantum force and local energy computations.
class FermionNInumerical : public WaveFunction
{
    public:
        // Constructor of the FermionNInumerical class. Initializes the parameters for the fermion system.
        //
        // Input:   double alpha - variational parameter alpha;
        //          int n_particles - number of particles in the system;
        //          double omega - angular frequency of the system.
        FermionNInumerical(double alpha, int n_particles, double omega);


        // Fills the state vector x with particle positions from the system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  VectorXvar - the state vector containing the particle positions.
        VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);


        // Computes the value of Psi1 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute Psi1.
        //
        // Output:  var - the value of Psi1 for the given particle.
        var psi1i(VectorXvar& x, int idx);


        // Computes the value of Psi2 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute Psi2.
        //
        // Output:  var - the value of Psi2 for the given particle.
        var psi2i(VectorXvar& x, int idx);


        // Computes the value of Psi3 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute Psi3.
        //
        // Output:  var - the value of Psi3 for the given particle.
        var psi3i(VectorXvar& x, int idx);


        // Computes the value of Psi4 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute Psi4.
        //
        // Output:  var - the value of Psi4 for the given particle.
        var psi4i(VectorXvar& x, int idx);


        // Computes the value of Psi5 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute Psi5.
        //
        // Output:  var - the value of Psi5 for the given particle.
        var psi5i(VectorXvar& x, int idx);


        // Computes the value of Psi6 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute Psi6.
        //
        // Output:  var - the value of Psi6 for the given particle.
        var psi6i(VectorXvar& x, int idx);


        // Computes the gradient of Psi1 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi1.
        //
        // Output:  VectorXvar - the gradient of Psi1 for the given particle.
        VectorXvar GradPsi1i(VectorXvar& x, int idx);


        // Computes the gradient of Psi2 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi2.
        //
        // Output:  VectorXvar - the gradient of Psi2 for the given particle.
        VectorXvar GradPsi2i(VectorXvar& x, int idx);


        // Computes the gradient of Psi3 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi3.
        //
        // Output:  VectorXvar - the gradient of Psi3 for the given particle.
        VectorXvar GradPsi3i(VectorXvar& x, int idx);


        // Computes the gradient of Psi4 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi4.
        //
        // Output:  VectorXvar - the gradient of Psi4 for the given particle.
        VectorXvar GradPsi4i(VectorXvar& x, int idx);


        // Computes the gradient of Psi5 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi5.
        //
        // Output:  VectorXvar - the gradient of Psi5 for the given particle.
        VectorXvar GradPsi5i(VectorXvar& x, int idx);


        // Computes the gradient of Psi6 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi6.
        //
        // Output:  VectorXvar - the gradient of Psi6 for the given particle.
        VectorXvar GradPsi6i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi1 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi1.
        //
        // Output:  var - the Laplacian of Psi1 for the given particle.
        var LaplPsi1i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi2 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi2.
        //
        // Output:  var - the Laplacian of Psi2 for the given particle.
        var LaplPsi2i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi3 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi3.
        //
        // Output:  var - the Laplacian of Psi3 for the given particle.
        var LaplPsi3i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi4 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi4.
        //
        // Output:  var - the Laplacian of Psi4 for the given particle.
        var LaplPsi4i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi5 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi5.
        //
        // Output:  var - the Laplacian of Psi5 for the given particle.
        var LaplPsi5i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi6 for the specified particle using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi6.
        //
        // Output:  var - the Laplacian of Psi6 for the given particle.
        var LaplPsi6i(VectorXvar& x, int idx);


        // Computes the Slater determinant or its derivatives (gradient or Laplacian) for the fermionic system using the given parameters.
        // The function computes the value of the Slater determinant, its gradient, or its Laplacian depending on the `der_order` argument.
        //
        // Input:   VectorXvar& x - state vector containing particle positions;
        //          int particles - number of particles in the system;
        //          int row_changed - row index of the particle that has been changed;
        //          int der_order - order of the derivatives to compute (0 - value, 1 - gradient, 2 - Laplacian);
        //          int grad_comp - gradient component to compute (if applicable, typically x or y component).
        //
        // Output:  double - the computed Slater determinant value or its derivative (gradient or Laplacian) depending on `der_order`.
        double SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp);
        

        // Computes the sum of Laplacians of Psi_T divided by Psi_T for all particles in the system using the fermionic ansatz.
        //
        // Input:   VectorXvar& x - state vector containing particle positions.
        //
        // Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles.
        double LaplPsiTOverPsiT(VectorXvar& x);


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


        // WaveFunction class function that computes the quantum force for a specific particle in the system
        // using the fermionic ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int n - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the quantum force for the specified particle.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
};