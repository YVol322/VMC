#pragma once

#include <memory>
#include <vector>
#include <cassert>

#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"


class System {
public:
    System(
            std::unique_ptr<class Hamiltonian> hamiltonian,
            std::unique_ptr<class WaveFunction> waveFunction,
            std::unique_ptr<class MonteCarlo> solver,
            std::vector<std::unique_ptr<class Particle>> particles);

    unsigned int runEquilibrationSteps(
            double stepLength,
            unsigned int numberOfEquilibrationSteps);

    std::unique_ptr<class Sampler> runMetropolisSteps(
            double stepLength,
            unsigned int numberOfMetropolisSteps);

    double computeLocalEnergy();
    double computerij(int i, int j);
    const std::vector<double>& getWaveFunctionParameters();

private:
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;

    std::unique_ptr<class Hamiltonian> m_hamiltonian;
    std::unique_ptr<class WaveFunction> m_waveFunction;
    std::unique_ptr<class MonteCarlo> m_solver;
    std::vector<std::unique_ptr<class Particle>> m_particles;
};