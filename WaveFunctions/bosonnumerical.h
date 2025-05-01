#pragma once

#include "wavefunction.h"

#include <cmath>
#include <iostream>

#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
using namespace autodiff;
using autodiff::var;
using autodiff::VectorXvar;
using autodiff::MatrixXvar;

class BosonNumerical : public WaveFunction
{
    public:
        // Constructor.
        BosonNumerical(double alpha, std::vector<double> beta, int mode, int n_particles);

        // BosonNumerical subclass specific functions.
        VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);
        int BetaIndex(int i, int j);
        var a_ij(int i, int j);
        var Jastrow(VectorXvar& x);
        var PadeJastrow(VectorXvar& x);
        var Phi(VectorXvar& x);
        VectorXvar dJdxi(VectorXvar& x);
        VectorXvar dPhidxi(VectorXvar& x);
        VectorXvar d2Phidxi2(VectorXvar& x);
        VectorXvar d2Jdxi2(VectorXvar& x);
        double LaplPhiOverPhi(VectorXvar& x);
        double LaplJOverJ(VectorXvar& x);
        double GradGrad(VectorXvar& x);

        // WaveFunction class functions.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);
};