#include <cmath>
//#include <armadillo>
#include <Eigen/Dense>
#include <iostream>


#include "fermionsjastrow2autodiff.h"

FermionsJastrow2Autodiff::FermionsJastrow2Autodiff(double alpha, double beta)
{
    assert(alpha >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(m_numberOfParameters);

    m_parameters.push_back(beta);
    m_parameters.push_back(alpha);
}

var FermionsJastrow2Autodiff::Psi1(std::vector<std::unique_ptr<class Particle>>& particles)
{
    var alpha = m_parameters.back();
    int n_particles = particles.size();
    
    var sum = 0;
    for(int i = 0; i < n_particles; i++)
    {
        sum += r_squared(particles, i);
    }

    sum *= -alpha;

    return exp(sum);
}

var FermionsJastrow2Autodiff::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    var beta12 = m_parameters.at(0);
    Particle& p1 = *(particles.at(0));
    Particle& p2 = *(particles.at(1));

    var x1 = p1.getPosition().at(0);
    var y1 = p1.getPosition().at(1);
    var x2 = p2.getPosition().at(0);
    var y2 = p2.getPosition().at(1);

    var dx = x1 - x2;
    var dy = y1 - y2;
    var r12 = sqrt(dx * dx + dy * dy);

    return exp(beta12 * r12);
}



double FermionsJastrow2Autodiff::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double psi1 = val(Psi1(particles));
    double J = val(Jastrow(particles));

    return psi1 * J;
}

double FermionsJastrow2Autodiff::GradPsi1GradJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles)
{
    var alpha = m_parameters.at(1);
    var beta = m_parameters.at(0);

    // Define autodiff vars
    var x1 = particles[0]->getPosition().at(0);
    var y1 = particles[0]->getPosition().at(1);
    var x2 = particles[1]->getPosition().at(0);
    var y2 = particles[1]->getPosition().at(1);

    // Construct psi1 from autodiff vars
    auto psi1_func = [&](var x1_, var y1_, var x2_, var y2_) {
        return exp(-alpha * (x1_ * x1_ + y1_ * y1_ + x2_ * x2_ + y2_ * y2_));
    };

    // Construct Jastrow from autodiff vars
    auto J_func = [&](var x1_, var y1_, var x2_, var y2_) {
        var dx = x1_ - x2_;
        var dy = y1_ - y2_;
        return exp(beta * sqrt(dx * dx + dy * dy));
    };

    var psi1 = psi1_func(x1, y1, x2, y2);
    var J = J_func(x1, y1, x2, y2);

    // Gradients of both
    auto [psi1_x1, psi1_y1, psi1_x2, psi1_y2] = derivativesx(psi1, wrt(x1, y1, x2, y2));
    auto [J_x1, J_y1, J_x2, J_y2] = derivativesx(J, wrt(x1, y1, x2, y2));

    var cross_term = 2.0 * (psi1_x1 * J_x1 + psi1_y1 * J_y1 + psi1_x2 * J_x2 + psi1_y2 * J_y2);

    return val(cross_term) / (val(psi1) * val(J));
}


double FermionsJastrow2Autodiff::LaplacianPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles)
{
    var alpha = m_parameters.at(1);

    auto psi1 = [&](var x1, var y1, var x2, var y2)
    {
        var alpha = m_parameters.at(1);
        var r2 = x1*x1 + y1*y1 + x2*x2 + y2*y2;
        return exp(-alpha * r2);
    };

    var x1 = particles[0]->getPosition().at(0);
    var y1 = particles[0]->getPosition().at(1);
    var x2 = particles[1]->getPosition().at(0);
    var y2 = particles[1]->getPosition().at(1);

    // now compute value and derivatives
    var psi_val = psi1(x1, y1, x2, y2);

    auto [dx1, dy1, dx2, dy2] = derivativesx(psi_val, wrt(x1, y1, x2, y2));
    auto [dxx1] = derivativesx(dx1, wrt(x1));
    auto [dyy1] = derivativesx(dy1, wrt(y1));
    auto [dxx2] = derivativesx(dx2, wrt(x2));
    auto [dyy2] = derivativesx(dy2, wrt(y2));

    double laplacian = val(dxx1 + dyy1 + dxx2 + dyy2);
    double psi_scalar = val(psi_val);

    return laplacian / psi_scalar;

}

double FermionsJastrow2Autodiff::LaplacianJOverJ(std::vector<std::unique_ptr<class Particle>>& particles)
{
    var beta = m_parameters.at(0);

    // Define autodiff variables
    var x1 = particles[0]->getPosition().at(0);
    var y1 = particles[0]->getPosition().at(1);
    var x2 = particles[1]->getPosition().at(0);
    var y2 = particles[1]->getPosition().at(1);

    // Define Jastrow directly in terms of autodiff variables
    auto JastrowFunc = [&](var x1_, var y1_, var x2_, var y2_) {
        var dx = x1_ - x2_;
        var dy = y1_ - y2_;
        var r12 = sqrt(dx * dx + dy * dy);
        return exp(beta * r12);
    };

    var J = JastrowFunc(x1, y1, x2, y2);

    // First derivatives
    auto [J_x1, J_y1, J_x2, J_y2] = derivativesx(J, wrt(x1, y1, x2, y2));

    // Second derivatives (Laplacian)
    auto [J_x1x1] = derivativesx(J_x1, wrt(x1));
    auto [J_y1y1] = derivativesx(J_y1, wrt(y1));
    auto [J_x2x2] = derivativesx(J_x2, wrt(x2));
    auto [J_y2y2] = derivativesx(J_y2, wrt(y2));

    double laplacian = val(J_x1x1) + val(J_y1y1) + val(J_x2x2) + val(J_y2y2);

    return laplacian / val(J);
}



double FermionsJastrow2Autodiff::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double laplacianPsi1 = LaplacianPsi1OverPsi1(particles);
    double laplacianJ = LaplacianJOverJ(particles);
    double crossTerm = GradPsi1GradJOverPsi(particles);

    double result = laplacianPsi1 + laplacianJ + crossTerm;

    return result;
}

double FermionsJastrow2Autodiff::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
{
    Particle particle_i = *(particles.at(part_index));
    int n_dimensions = particle_i.getNumberOfDimensions();
    double coordinate;
    double r2 = 0;

    for(int i = 0; i < n_dimensions; i++)
    {
        coordinate = particle_i.getPosition().at(i);
        r2 += coordinate * coordinate;
    }

    return r2;
}

double FermionsJastrow2Autodiff::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
{
    Particle& pi = *(particles.at(i));
    Particle& pj = *(particles.at(j));

    double x_i = pi.getPosition().at(0);
    double y_i = pi.getPosition().at(1);
    double x_j = pj.getPosition().at(0);
    double y_j = pj.getPosition().at(1);

    double r_ij = sqrt((x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j));

    return r_ij;
}