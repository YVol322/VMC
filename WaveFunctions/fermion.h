#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

#include "wavefunction.h"

class Fermion : public WaveFunction
{
    public:
        Fermion(double alpha, std::vector<double> beta, int mode, int n_particles);
    
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);
        double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);
    
    
        double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        std::vector<double> GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        std::vector<double> GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        std::vector<double> GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        std::vector<double> GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        std::vector<double> GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        std::vector<double> GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    
        double SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp);
    
        double LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles);
        double GradiPsi1GradiJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
    
        double PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
};