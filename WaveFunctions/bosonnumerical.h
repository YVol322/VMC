#pragma once

#include <cmath>                            // Include the C++ math library for mathematical functions.
#include <iostream>                         // Include the C++ input-output stream library.

#include "wavefunction.h"                   // Include "wavefunction" header file with declarations.

#pragma once

#include <autodiff/reverse/var.hpp>         // Include autodiff reverse-mode library for automatic differentiation.
#include <autodiff/reverse/var/eigen.hpp>   // Include Eigen support for autodiff reverse-mode variables.
using namespace autodiff;
using autodiff::var;
using autodiff::VectorXvar;
using autodiff::MatrixXvar;


// Declaration of the BosonNumerical class. This class is a subclass of the WaveFunction class.
// Its purpose is to implement the evaluate, computeDoubleDerivative, and quantumForce functions 
// using automatic differentiation for the Bosonic ansatz.
// All other functions are designed to simplify these three core functions.
class BosonNumerical : public WaveFunction
{
    public:
        // Constructor of the BosonNumerical class. Initializes the parameters for the boson system.
        //
        // Input:   double alpha - variational parameter alpha;
        //          int n_particles - number of particles in the system;
        //          double omega - angular frequency of the system.
        BosonNumerical(double alpha, int n_particles, double omega);


        // BosonNumerical subclass specific function that fills the state vector x with the particle positions.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  VectorXvar - the vector of particle positions (state vector x).
        VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);


        // BosonNumerical subclass specific function that computes the trial wave function (Psi_T) using the state vector x.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  var - the value of the trial wave function (Psi_T) for the given particle configuration.
        var PsiT(VectorXvar& x);


        // BosonNumerical subclass specific function that computes the Laplacian of Psi_T over Psi_T using the state vector x
        // and automatic differentiation.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions.
        //
        // Output:  double - the computed Laplacian of Psi_T divided by Psi_T for all particles in the system.
        double LaplPsiTOverPsiT(VectorXvar& x);


        // WaveFunction class function that evaluates the wave function for the given particle configuration 
        // using the bosonic ansatz.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the value of the wave function for the given particle configuration using the bosonic ansatz.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
        // for the bosonic ansatz, which is required for calculating the local energy using automatic differentiation.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);


        // WaveFunction class function that computes the quantum force for a specific particle in the system
        // using the bosonic ansatz and automatic differentiation.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int n - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the quantum force for the specified particle.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
};