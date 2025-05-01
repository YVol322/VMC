#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

#include "wavefunction.h"

class Fermion : public WaveFunction
{
    public:
        // Constructor.
        Fermion(double alpha, std::vector<double> beta, int mode, int n_particles);
    
        // Fermion subclass specific functions.
        int BetaIndex(int i, int j);
        double a_ij(int i, int j);
        double Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        double PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles);
        double Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp);
        std::vector<double> GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        std::vector<double> GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);
        double LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx);
        double LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles);
        double GradiPsi1GradiJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx);

        // WaveFunction class functions.
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx);
        double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);
};