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
        BosonNumerical(double alpha, int n_particles);

        // BosonNumerical subclass specific functions.
        VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);
        var PsiT(VectorXvar& x);
        double LaplPsiTOverPsiT(VectorXvar& x);

        // WaveFunction class functions.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
};