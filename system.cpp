#include <iostream>

#include "system.h"


System::System(
        std::unique_ptr<class Hamiltonian> hamiltonian,
        std::unique_ptr<class WaveFunction> waveFunction,
        std::unique_ptr<class MonteCarlo> solver,
        std::vector<std::unique_ptr<class Particle>> particles)
{
    m_numberOfParticles = particles.size();;
    m_numberOfDimensions = particles[0]->getNumberOfDimensions();
    m_hamiltonian = std::move(hamiltonian);
    m_waveFunction = std::move(waveFunction);
    m_solver = std::move(solver);
    m_particles = std::move(particles);
}


unsigned int System::runEquilibrationSteps(
        double stepLength,
        unsigned int numberOfEquilibrationSteps)
{
    unsigned int acceptedSteps = 0;

    for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++)
    {
        acceptedSteps += m_solver->step(stepLength, *m_waveFunction, m_particles);
    }

    return acceptedSteps;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
        double stepLength,
        unsigned int numberOfMetropolisSteps)
{
    auto sampler = std::make_unique<Sampler>(
            m_numberOfParticles,
            m_numberOfDimensions,
            stepLength,
            numberOfMetropolisSteps);

    for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
        /* Call solver method to do a single Monte-Carlo step.
         */
        bool acceptedStep = m_solver->step(stepLength, *m_waveFunction, m_particles);

        /* Here you should sample the energy (and maybe other things) using the
         * sampler instance of the Sampler class.
         */
        sampler->sample(acceptedStep, this);
    }

    sampler->computeAverages();

    return sampler;
}

double System::computeLocalEnergy()
{
    // Helper function
    return m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles);
}

double System::computerij(int i, int j)
{
    // Helper function
    return m_waveFunction -> r_ij(m_particles, i, j);
}

double System::computer2()
{
    double r2 = 0;
    for(int i = 0; i < m_numberOfParticles; i++)
    {
        r2 += m_waveFunction -> r_squared(m_particles, i);
    }

    return r2;
}

const std::vector<double>& System::getWaveFunctionParameters()
{
    // Helper function
    return m_waveFunction->getParameters();
}