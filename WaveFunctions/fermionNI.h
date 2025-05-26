#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <iostream>

#include "wavefunction.h"

class FermionNI : public WaveFunction
{
    public:
        // Constructor.
        FermionNI(double alpha, int n_particles);
    
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
        double LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles);
        double SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp);

        // WaveFunction class functions.
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n);
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
};