#pragma once

#include "wavefunction.h"

#include <autodiff/reverse/var.hpp>
using namespace autodiff;
using autodiff::var;

class FermionsJastrow2Autodiff : public WaveFunction {
public:
    FermionsJastrow2Autodiff(double alpha, double beta);
    
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);
    double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);

    var Psi1(std::vector<std::unique_ptr<class Particle>>& particles);
    var Jastrow(std::vector<std::unique_ptr<class Particle>>& particles);

    double GradPsi1GradJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles);
    double LaplacianPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles);
    double LaplacianJOverJ(std::vector<std::unique_ptr<class Particle>>& particles);

};