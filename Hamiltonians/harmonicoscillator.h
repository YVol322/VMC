#pragma once

#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian
{
    public:
        HarmonicOscillator(double omega, bool Coulomb);
        double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        );

    private:
        double m_omega;
        bool m_Coulomb;
};