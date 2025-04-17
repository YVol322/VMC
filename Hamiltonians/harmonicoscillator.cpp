#include <iostream>

#include "harmonicoscillator.h"

HarmonicOscillator::HarmonicOscillator(double omega)
{
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        )
{
    int n_particles = particles.size();
    double r2 = 0.0;

    for (int i = 0; i < n_particles; i++) {
        r2 += waveFunction.r_squared(particles, i);
    }

    double potentialEnergy = 0.5 * r2 * m_omega * m_omega;
    double kineticEnergy = -0.5 * waveFunction.computeDoubleDerivative(particles);

    // Coulomb repulsion
    double coulombEnergy = 0.0;
    for (int i = 0; i < n_particles - 1; ++i) {
        for (int j = i + 1; j < n_particles; ++j) {
            double rij = computerij(waveFunction, particles, i, j);
            coulombEnergy += 1.0 / rij;
        }
    }

    return kineticEnergy + potentialEnergy + coulombEnergy;
}


double HarmonicOscillator::computerij(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles,
            int i, int j
        )
{
    Particle pi = *(particles.at(i));
    Particle pj = *(particles.at(j));

    double x_i = pi.getPosition().at(0);
    double y_i = pi.getPosition().at(1);
    double x_j = pj.getPosition().at(0);
    double y_j = pj.getPosition().at(1);

    double r_ij = sqrt((x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j));

    return r_ij;
}