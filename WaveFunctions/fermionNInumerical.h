#pragma once

#include "wavefunction.h"

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
using namespace autodiff;
using autodiff::var;
using autodiff::VectorXvar;
using autodiff::MatrixXvar;


class FermionNInumerical : public WaveFunction
{
    public:
        // Constructor.
        FermionNInumerical(double alpha, int n_particles, double omega);
    
        VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);
        var psi1i(VectorXvar& x, int idx);
        var psi2i(VectorXvar& x, int idx);
        var psi3i(VectorXvar& x, int idx);
        var psi4i(VectorXvar& x, int idx);
        var psi5i(VectorXvar& x, int idx);
        var psi6i(VectorXvar& x, int idx);
        VectorXvar GradPsi1i(VectorXvar& x, int idx);
        VectorXvar GradPsi2i(VectorXvar& x, int idx);
        VectorXvar GradPsi3i(VectorXvar& x, int idx);
        VectorXvar GradPsi4i(VectorXvar& x, int idx);
        VectorXvar GradPsi5i(VectorXvar& x, int idx);
        VectorXvar GradPsi6i(VectorXvar& x, int idx);
        var LaplPsi1i(VectorXvar& x, int idx);
        var LaplPsi2i(VectorXvar& x, int idx);
        var LaplPsi3i(VectorXvar& x, int idx);
        var LaplPsi4i(VectorXvar& x, int idx);
        var LaplPsi5i(VectorXvar& x, int idx);
        var LaplPsi6i(VectorXvar& x, int idx);
        double SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp);
        double LaplPsiTOverPsiT(VectorXvar& x);

        // WaveFunction class functions.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
};