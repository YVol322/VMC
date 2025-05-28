#include <iostream>

#include "harmonicoscillator.h"

HarmonicOscillator::HarmonicOscillator(double omega, bool Coulomb)
{
    assert(omega > 0);
    m_omega  = omega;
    m_Coulomb = Coulomb;
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

    double localEnergy = potentialEnergy + kineticEnergy;

    if(m_Coulomb)
    {
        double rij;
        double coulombEnergy = 0.0;
        for (int i = 0; i < n_particles; ++i) {
            for (int j = i + 1; j < n_particles; ++j) {
                rij = waveFunction.r_ij(particles, i, j);
                coulombEnergy += 1.0 / rij;
            }
        }
        localEnergy += coulombEnergy;
    }

    return localEnergy;
}