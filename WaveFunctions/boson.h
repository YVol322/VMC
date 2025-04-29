#pragma once

#include "wavefunction.h"

class Boson : public WaveFunction
{
    public:
        Boson(double alpha, std::vector<double> beta, int mode, int n_particles);
        double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);

        double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        double PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx);
        std::vector<double> GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx);
        double LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx);
        double LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx);
        double LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles);
        double GradPsi1GradJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles);
};