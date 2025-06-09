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


class FermionNumerical : public WaveFunction
{
    public:
        // Constructor of the FermionNumerical class. Initializes the parameters for the fermion system.
        //
        // Input:   double alpha - variational parameter alpha;
        //          std::vector<double> beta - vector of variational parameters for Jastrow or Pade-Jastrow;
        //          int mode - 0 for Jastrow ansatz, 1 for Pade-Jastrow ansatz;
        //          int n_particles - number of particles in the system;
        //          double omega - angular frequency of the system.
        FermionNumerical(double alpha, std::vector<double> beta, int mode, int n_particles, double omega);


        // Fills the state vector x with the particle positions for automatic differentiation.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  VectorXvar - state vector containing the particle positions.
        VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);


        // Computes the Jastrow factor for the fermionic system using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  var - the computed Jastrow factor for the given particle configuration.
        var Jastrow(VectorXvar& x);


        // Computes the Pade-Jastrow factor for the fermionic system using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  var - the computed Pade-Jastrow factor for the given particle configuration.
        var PadeJastrow(VectorXvar& x);


        // Computes the value of Psi1 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute Psi1.
        //
        // Output:  var - the value of Psi1 for the given particle.
        var psi1i(VectorXvar& x, int idx);


        // Computes the value of Psi2 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute Psi2.
        //
        // Output:  var - the value of Psi2 for the given particle.
        var psi2i(VectorXvar& x, int idx);


        // Computes the value of Psi3 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute Psi3.
        //
        // Output:  var - the value of Psi3 for the given particle.
        var psi3i(VectorXvar& x, int idx);


        // Computes the value of Psi4 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute Psi4.
        //
        // Output:  var - the value of Psi4 for the given particle.
        var psi4i(VectorXvar& x, int idx);


        // Computes the value of Psi5 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute Psi5.
        //
        // Output:  var - the value of Psi5 for the given particle.
        var psi5i(VectorXvar& x, int idx);


        // Computes the value of Psi6 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute Psi6.
        //
        // Output:  var - the value of Psi6 for the given particle.
        var psi6i(VectorXvar& x, int idx);


        // Computes the gradient of Psi1 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi1.
        //
        // Output:  VectorXvar - the gradient of Psi1 for the given particle.
        VectorXvar GradPsi1i(VectorXvar& x, int idx);


        // Computes the gradient of Psi2 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi2.
        //
        // Output:  VectorXvar - the gradient of Psi2 for the given particle.
        VectorXvar GradPsi2i(VectorXvar& x, int idx);


        // Computes the gradient of Psi3 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi3.
        //
        // Output:  VectorXvar - the gradient of Psi3 for the given particle.
        VectorXvar GradPsi3i(VectorXvar& x, int idx);


        // Computes the gradient of Psi4 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi4.
        //
        // Output:  VectorXvar - the gradient of Psi4 for the given particle.
        VectorXvar GradPsi4i(VectorXvar& x, int idx);


        // Computes the gradient of Psi5 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi5.
        //
        // Output:  VectorXvar - the gradient of Psi5 for the given particle.
        VectorXvar GradPsi5i(VectorXvar& x, int idx);


        // Computes the gradient of Psi6 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of Psi6.
        //
        // Output:  VectorXvar - the gradient of Psi6 for the given particle.
        VectorXvar GradPsi6i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi1 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi1.
        //
        // Output:  var - the Laplacian of Psi1 for the given particle.
        var LaplPsi1i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi2 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi2.
        //
        // Output:  var - the Laplacian of Psi2 for the given particle.
        var LaplPsi2i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi3 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi3.
        //
        // Output:  var - the Laplacian of Psi3 for the given particle.
        var LaplPsi3i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi4 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi4.
        //
        // Output:  var - the Laplacian of Psi4 for the given particle.
        var LaplPsi4i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi5 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi5.
        //
        // Output:  var - the Laplacian of Psi5 for the given particle.
        var LaplPsi5i(VectorXvar& x, int idx);


        // Computes the Laplacian of Psi6 for the specified particle using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of Psi6.
        //
        // Output:  var - the Laplacian of Psi6 for the given particle.
        var LaplPsi6i(VectorXvar& x, int idx);


        // Computes the Slater determinant or its gradient component or Laplacian depending on the derivative order.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int particles - the particle set (spin-up or spin-down);
        //          int row_changed - row index of the particle that was changed;
        //          int der_order - order of the derivatives to compute (0 for Psi, 1 for gradient, 2 for Laplacian);
        //          int grad_comp - component of the gradient to compute (if applicable).
        //
        // Output:  double - the computed value based on the requested derivative order.
        double SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp);


        // Computes the gradient of the Jastrow factor for a given particle.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of the Jastrow factor.
        //
        // Output:  VectorXvar - the gradient of the Jastrow factor for the given particle.
        VectorXvar GradiJastrow(VectorXvar& x, int idx);


        // Computes the gradient of the Pade-Jastrow factor for a given particle.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the gradient of the Pade-Jastrow factor.
        //
        // Output:  VectorXvar - the gradient of the Pade-Jastrow factor for the given particle.
        VectorXvar GradiPadeJastrow(VectorXvar& x, int idx);


        // Computes the Laplacian of the Jastrow factor for a given particle.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of the Jastrow factor.
        //
        // Output:  var - the Laplacian of the Jastrow factor for the given particle.
        var LapliJastrow(VectorXvar& x, int idx);


        // Computes the Laplacian of the Pade-Jastrow factor for a given particle.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int idx - index of the particle for which to compute the Laplacian of the Pade-Jastrow factor.
        //
        // Output:  var - the Laplacian of the Pade-Jastrow factor for the given particle.
        var LapliPadeJastrow(VectorXvar& x, int idx);


        // Computes the Laplacian of PsiT divided by PsiT for all particles in the system using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  double - the computed sum of Laplacians of PsiT divided by PsiT for all particles.
        double LaplPsiTOverPsiT(VectorXvar& x);


        // Computes the Laplacian of the Jastrow factor divided by the Jastrow factor for all particles in the system using automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  double - the computed sum of Laplacians of Jastrow divided by Jastrow for all particles.
        double LaplJOverJ(VectorXvar& x);


        // Computes the gradient of PsiT multiplied by the gradient of the Jastrow factor divided by PsiT for all particles in the system.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  double - the computed sum of the gradient of PsiT multiplied by the gradient of Jastrow over PsiT for all particles.
        double GradPsiTGradJOverPsiT(VectorXvar& x);



        // Computes the value of the wave function for the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed value of the wave function for the given particle configuration.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);


        // Computes the sum of second derivatives (Laplace operator) of the wave function for the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);


        // Computes the quantum force for a specific particle in the fermionic system.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int i - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the quantum force for the specified particle.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int i);
};