#include <memory>
#include <cmath>
#include <cassert>

#include "simplegaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

SimpleGaussian::SimpleGaussian(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
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

double SimpleGaussian::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    double alpha = m_parameters.back();
    int n_particles = particles.size();
    double argument = 0;
    double wavefunction;

    for(int i = 0; i < n_particles; i++)
    {
        argument += r_squared(particles, i);
    }

    argument *= -0.5 * alpha * alpha;
    wavefunction = exp(argument);

    return wavefunction;
}

double SimpleGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    int n_particles = particles.size();
    double summ_of_squares = 0;
    double laplacian;

    for(int i = 0; i < n_particles; i++)
    {
        summ_of_squares += r_squared(particles, i);
    }

    int n_dims = (particles.back()) -> getNumberOfDimensions();
    laplacian = - n_dims * n_particles * alpha * alpha + alpha * alpha * alpha * alpha * summ_of_squares;

    return laplacian;
}