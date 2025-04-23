#pragma once

#include "wavefunction.h"

class FermionsJastrow : public WaveFunction {
public:
    FermionsJastrow(double alpha, std::vector<double> beta);

    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);
    double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);

    // Новые функции будут добавляться позже
    double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);
    std::vector<double> GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    std::vector<double> GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    std::vector<double> GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    std::vector<double> GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);

    double SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp);

    double LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles);
    double GradiPsi1GradiJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);



};
