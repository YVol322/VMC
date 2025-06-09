#include <iostream>                 // Include the C++ input-output stream library.

#include "harmonicoscillator.h"     // Include "harmonicoscillator" header file with declarations.


// Constructor of the HarmonicOscillator class. Initializes the angular frequency (omega) and the Coulomb flag.
//
// Input:   double omega - the angular frequency of the harmonic oscillator;
//          bool Coulomb - indicates whether Coulomb interaction is included.
HarmonicOscillator::HarmonicOscillator(double omega, bool Coulomb)
{
    assert(omega > 0);
    m_omega  = omega;
    m_Coulomb = Coulomb;
}


// Function that computes the local energy for a system of particles in the harmonic oscillator potential.
//
// Input:   class WaveFunction& waveFunction - the wavefunction of the system;
//          std::vector<std::unique_ptr<class Particle>>& particles - the list of particles in the system.
//
// Output:  double - the computed local energy.
double HarmonicOscillator::computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        )
{
    int n_particles = particles.size();

    double r2 = 0.0;    // Initialize the variable r2 to store the sum of the squared coordinates of all particles.

    // Calculate the sum of squared coordinates for all particles and store it in r2.
    for (int i = 0; i < n_particles; i++)
    {
        r2 += waveFunction.r_squared(particles, i);
    }

    // Harmonic oscillator potential term in the local energy.
    double potentialEnergy = 0.5 * r2 * m_omega * m_omega;

    // Kinetic energy term in the local energy.
    double kineticEnergy = -0.5 * waveFunction.computeDoubleDerivative(particles);


    double localEnergy = potentialEnergy + kineticEnergy;

    // If m_Coulomb is true, add the Coulomb term to the local energy.
    if(m_Coulomb)
    {
        double rij;
        double coulombEnergy = 0.0;
        for (int i = 0; i < n_particles; i++)
        {
            for (int j = i + 1; j < n_particles; j++)
            {
                rij = waveFunction.r_ij(particles, i, j);
                coulombEnergy += 1.0 / rij;
            }
        }
        localEnergy += coulombEnergy;
    }

    return localEnergy;
}