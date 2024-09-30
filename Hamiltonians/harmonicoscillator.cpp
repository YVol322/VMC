#include<memory>
#include <cassert>
#include <iostream>

#include "harmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

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
    double r2 = 0;

    //#pragma omp parallel for
    for(int i = 0; i < n_particles; i++)
    {
        r2 += waveFunction.r_squared(particles, i);
    }


    double potentialEnergy = 0.5 * r2 * m_omega * m_omega ;

    // Computing kinetic energy (without exponent).
    double kineticEnergy   = -0.5 * waveFunction.computeDoubleDerivative(particles);

    // Computing local energy (based on the analytical expression).
    // Note! I computed kinetic & potential energy without exponent (wave fucntion),
    // because it cancles after division. This was made to fasten program.
    return kineticEnergy + potentialEnergy;
}
