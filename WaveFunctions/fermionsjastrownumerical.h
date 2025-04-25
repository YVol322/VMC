#pragma once

#include "wavefunction.h"

#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
using namespace autodiff;
using autodiff::var;
using autodiff::VectorXvar;
using autodiff::MatrixXvar;

class FermionsJastrowNumerical : public WaveFunction {
public:
    FermionsJastrowNumerical(double alpha, std::vector<double> beta);
    
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
    double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);
    double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);

    VectorXvar fill_x(std::vector<std::unique_ptr<class Particle>>& particles);

    var psi1i(VectorXvar& x, int idx);
    var psi2i(VectorXvar& x, int idx);
    var psi3i(VectorXvar& x, int idx);
    var psi4i(VectorXvar& x, int idx);
    var psi5i(VectorXvar& x, int idx);
    var psi6i(VectorXvar& x, int idx);
    VectorXvar GradPsi1i(VectorXvar& x, int idx);
    VectorXvar GradPsi2i(VectorXvar& x, int idx);
    VectorXvar GradPsi3i(VectorXvar& x, int idx);
    VectorXvar GradPsi4i(VectorXvar& x, int idx);
    VectorXvar GradPsi5i(VectorXvar& x, int idx);
    VectorXvar GradPsi6i(VectorXvar& x, int idx);
    var LaplPsi1i(VectorXvar& x, int idx);
    var LaplPsi2i(VectorXvar& x, int idx);
    var LaplPsi3i(VectorXvar& x, int idx);
    var LaplPsi4i(VectorXvar& x, int idx);
    var LaplPsi5i(VectorXvar& x, int idx);
    var LaplPsi6i(VectorXvar& x, int idx);
    var Jastrow(const VectorXvar& x);
    VectorXvar GradiJastrow(VectorXvar& x, int idx);
    var LapliJastrow(VectorXvar& x, int idx);

    double SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp);

    double GradPsi1GradJOverPsi(VectorXvar& x);
    double LaplacianPsi1OverPsi1(VectorXvar& x);
    double LaplacianJOverJ(VectorXvar& x);

};