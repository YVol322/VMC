#pragma once

#include <cassert>
#include <memory>
#include <vector>

#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

class Hamiltonian
{
    public:
        virtual ~Hamiltonian() = default;
        virtual double computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        ) = 0;
};