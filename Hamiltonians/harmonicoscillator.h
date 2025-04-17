#pragma once

#include "hamiltonian.h"

class HarmonicOscillator : public Hamiltonian
{
    public:
        HarmonicOscillator(double omega);
        double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        );

        double computerij(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles,
            int i, int j
        );

    private:
        double m_omega;
};