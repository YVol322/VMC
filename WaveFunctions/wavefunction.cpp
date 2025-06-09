#include "wavefunction.h"    // Include "wavefunction" header file with declarations.


// Computes the sum of squared coordinates for the particle specified by part_idx.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int part_idx - index of the particle for which to compute r^2.
//
// Output:  double - the sum of squared coordinates for the specified particle.
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


// Computes the relative distance between particles i and j.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int i - index of the first particle;
//          int j - index of the second particle.
//
// Output:  double - the relative distance between particles i and j.
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


// Transforms a 2D index (i, j) into a 1D index k.
//
// Input:   int i - index of the first particle;
//          int j - index of the second particle.
//
// Output:  int - the corresponding 1D index.
int WaveFunction::BetaIndex(int i, int j)
{
    int n_particles = m_particles;

    int ii = std::min(i, j);
    int jj = std::max(i, j);

    return ii * (2 * n_particles - ii - 1) / 2 + (jj - ii - 1);
}


// Computes the constant a for particles i and j.
//
// Input:   int i - index of the first particle;
//          int j - index of the second particle.
//
// Output:  double - the constant a for particles i and j.
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


// Computes the Green's function ratio between two particle configurations.
//
// Input:   std::vector<double>& Rnew - new configuration of particle positions;
//          std::vector<double>& Rold - old configuration of particle positions;
//          double dt - time step used for the update;
//          std::vector<double>& Fold - forces for the old configuration;
//          std::vector<double>& Fnew - forces for the new configuration.
//
// Output:  double - the computed Green's function ratio.
double WaveFunction::GreensFunctionRatio(std::vector<double>& Rnew, std::vector<double>& Rold, 
                                         double dt, std::vector<double>& Fold, std::vector<double>& Fnew)
{
    const double D = 0.5;   // Constant D in a.u.
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