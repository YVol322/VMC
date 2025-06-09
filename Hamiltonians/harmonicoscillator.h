#pragma once

#include "hamiltonian.h"    // Include "hamiltonian" header file with declarations.


// Declaration of the HarmonicOscillator class. This class is a subclass of the Hamiltonian class,
// designed to compute the local energy for a system in a harmonic oscillator potential,
// optionally including Coulomb interactions.
//
//
// Private variables:    double m_omega - the angular frequency of the harmonic oscillator;
//                       bool m_Coulomb - flag indicating whether Coulomb interaction is included.
//
//
// Functions that objects of this class can use:    HarmonicOscillator(...);
//                                                  computeLocalEnergy(...);
class HarmonicOscillator : public Hamiltonian
{
    public:
        // Constructor of the HarmonicOscillator class. Initializes the angular frequency (omega) and the Coulomb flag.
        //
        // Input:   double omega - the angular frequency of the harmonic oscillator;
        //          bool Coulomb - indicates whether Coulomb interaction is included.
        HarmonicOscillator(double omega, bool Coulomb);


        // Function that computes the local energy for a system of particles in the harmonic oscillator potential.
        //
        // Input:   class WaveFunction& waveFunction - the wavefunction of the system;
        //          std::vector<std::unique_ptr<class Particle>>& particles - the list of particles in the system.
        //
        // Output:  double - the computed local energy.
        double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        );

    private:
        double m_omega;   // Angular frequency of the harmonic oscillator.
        bool m_Coulomb;   // Flag indicating whether Coulomb interaction is included.
};