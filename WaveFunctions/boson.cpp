#include "boson.h"


Boson::Boson(double alpha, int n_particles)
{
    assert(alpha >= 0);

    m_numberOfParameters = 1;

    m_parameters.reserve(m_numberOfParameters);
    m_parameters.push_back(alpha);

    m_particles = n_particles;
}


double Boson::PsiT(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    double argument = 0;

    for(int i = 0; i < m_particles; i++)
    {
        argument += r_squared(particles, i);
    }

    return exp(-alpha * argument);
}


double Boson::LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    double argument = 0;

    Particle& particle_i = *(particles[0]);
    int n_dimensions = particle_i.getNumberOfDimensions();

    for(int i = 0; i < m_particles; i++)
    {
        argument += r_squared(particles, i);
    }
    

    return (-2 * n_dimensions * m_particles * alpha + 4 * alpha * alpha * argument);
}


double Boson::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    return PsiT(particles);
}


double Boson::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    return LaplPsiTOverPsiT(particles);
}

std::vector<double> Boson::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    double alpha = m_parameters.back();
    double d = particles[n] -> getNumberOfDimensions();
    std::vector<double> r = particles[n] -> getPosition();
    std::vector<double> F(d);
    
    for(int i = 0; i < d; i++)
    {
        F[i] = -4 * alpha * r[i];
    }

    return F;
}