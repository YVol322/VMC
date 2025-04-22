#include <cmath>
//#include <armadillo>
#include <Eigen/Dense>
#include <iostream>


#include "fermionsjastrow2.h"

FermionsJastrow2::FermionsJastrow2(double alpha, double beta)
{
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(m_numberOfParameters);

    m_parameters.push_back(beta);
    m_parameters.push_back(alpha);
}


double FermionsJastrow2::Psi1(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    int n_particles = particles.size();
    
    double sum = 0;
    for(int i = 0; i < n_particles; i++)
    {
        sum += r_squared(particles, i);
    }

    sum *= -alpha;

    return exp(sum);
}


double FermionsJastrow2::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double beta12 = m_parameters.at(0);
    double r12 = r_ij(particles, 0, 1);

    return exp(beta12 * r12);
}


double FermionsJastrow2::Psi1LaplasianOverPsi(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    int n_particles = particles.size();
    
    double sum = 0;
    for(int i = 0; i < n_particles; i++)
    {
        sum += r_squared(particles, i);
    }

    double laplacian_psi1 = -8 * alpha + 4 * alpha * alpha * sum;

    return laplacian_psi1;
}


double FermionsJastrow2::GradPsiGradJast(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    double beta = m_parameters.at(0);
    double r12 = r_ij(particles, 0, 1);

    return -4 * alpha * beta * r12;
}


double FermionsJastrow2::laplacianJOverJ(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double beta12 = m_parameters.at(0);
    double r12 = r_ij(particles, 0, 1);

    return 2 * beta12 * (beta12 + 1/r12);
}

double FermionsJastrow2::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
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

double FermionsJastrow2::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
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

double FermionsJastrow2::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double psi1 = Psi1(particles);
    double J = Jastrow(particles);

    return psi1 * J;
}

double FermionsJastrow2::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double laplacian = 0;

    laplacian += Psi1LaplasianOverPsi(particles);
    laplacian += laplacianJOverJ(particles);
    laplacian += GradPsiGradJast(particles);

    return laplacian;
}