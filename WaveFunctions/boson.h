#pragma once

#include <cmath>
#include <iostream>

#include "wavefunction.h"

class Boson : public WaveFunction
{
    public:
        // Constructor.
        Boson(double alpha, int n_particles, double omega);

        // Boson subclass specific functions.
        double PsiT(std::vector<std::unique_ptr<class Particle>>& particles);
        double LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles);

        // WaveFunction class functions.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
};