#pragma once

#include <cmath>
#include <iostream>

#include "wavefunction.h"

class Boson : public WaveFunction
{
    public:
        // Constructor.
        Boson(double alpha, std::vector<double> beta, int mode, int n_particles);

        // Boson subclass specific functions.
        int BetaIndex(int i, int j);
        double a_ij(int i, int j);
        double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        double PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        double Phi(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        std::vector<double> GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        double LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        double LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        double LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles);
        double GradPsi1GradJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles);

        // WaveFunction class functions.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);
};