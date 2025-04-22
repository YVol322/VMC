#pragma once

#include "wavefunction.h"

class FermionsJastrow2 : public WaveFunction {
public:
    FermionsJastrow2(double alpha, double beta);

    double Psi1(std::vector<std::unique_ptr<class Particle>>& particles);
    double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);
    double Psi1LaplasianOverPsi(std::vector<std::unique_ptr<class Particle>>& particles);
    double GradPsiGradJast(std::vector<std::unique_ptr<class Particle>>& particles);
    double laplacianJOverJ(std::vector<std::unique_ptr<class Particle>>& particles);

    double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index) override;
    double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j) override;
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles) override;
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) override;

};
