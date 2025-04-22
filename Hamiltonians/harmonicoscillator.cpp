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
    double r12 = waveFunction.r_ij(particles, 0, 1);

    return kineticEnergy + potentialEnergy + 1.0 /r12 ;
}