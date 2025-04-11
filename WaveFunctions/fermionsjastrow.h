#pragma once

#include "wavefunction.h"

class FermionsJastrow : public WaveFunction {
public:
    FermionsJastrow(double alpha, std::vector<double> beta);

    // Реализуем наследуемые функции
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles) override;
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) override;
    double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index) override;

    // Новые функции будут добавляться позже
    double Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    std::vector<double> Psi1Der(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    std::vector<double> Psi2Der(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    std::vector<double> Psi3Der(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi1DoubleDer(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi2DoubleDer(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    double Psi3DoubleDer(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);

    double SD(std::vector<std::unique_ptr<class Particle>>& particles, double row_indx, double parts, int der_order, int coord);

    double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);

    std::vector<double> gradpsi(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
    std::vector<double> gradlogjast(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
    double gradpsijast(std::vector<std::unique_ptr<class Particle>>& particles);
    double lapllogjast(std::vector<std::unique_ptr<class Particle>>& particles);
    double jastrow(std::vector<std::unique_ptr<class Particle>>& particles);

};
