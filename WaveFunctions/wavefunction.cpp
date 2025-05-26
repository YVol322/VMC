#include "wavefunction.h"


double WaveFunction::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    Particle& particle_i = *(particles[part_idx]);
    int n_dimensions = particle_i.getNumberOfDimensions();
    double coordinate;
    double r2 = 0;

    for(int i = 0; i < n_dimensions; i++)
    {
        coordinate = particle_i.getPosition()[i];
        r2 += coordinate * coordinate;
    }

    return r2;
}


double WaveFunction::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
{
    Particle& pi = *(particles[i]);
    Particle& pj = *(particles[j]);

    int n_dimensions = pi.getNumberOfDimensions();

    double rij = 0;
    for(int d = 0; d < n_dimensions; d++)
    {
        rij += (pi.getPosition()[d] - pj.getPosition()[d]) * (pi.getPosition()[d] - pj.getPosition()[d]);
    }
    
    return sqrt(rij);
}


int WaveFunction::BetaIndex(int i, int j)
{
    int n_particles = m_particles;

    int ii = std::min(i, j);
    int jj = std::max(i, j);

    return ii * (2 * n_particles - ii - 1) / 2 + (jj - ii - 1);
}


double WaveFunction::a_ij(int i, int j)
{
    int n_particles = m_particles;

    double aij;
    
    if ((i < n_particles / 2 && j >= n_particles / 2) || (i >= n_particles / 2 && j < n_particles / 2))
    {
        aij = 1;
    }
    else
    {
        aij = 1.0 / 3.0;
    }

    return aij;
}


double WaveFunction::GreensFunctionRatio(std::vector<double>& Rnew, std::vector<double>& Rold, 
double dt, std::vector<double>& Fold, std::vector<double>& Fnew)
{
    const double D = 0.5;
    double sum = 0.0;

    int ndim = Rold.size();
    for (int j = 0; j < ndim; ++j)
    {
        double num = Rnew[j] - Rold[j] - D * Fold[j] * dt;
        double den = Rold[j] - Rnew[j] - D * Fnew[j] * dt;
        sum += (num * num) - (den * den);
    }

    return exp(sum / (4.0 * D * dt));
}