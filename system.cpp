#include <iostream>     // Include the C++ input-output stream library.

#include "system.h"     // Include "system" header file with declarations.


// Constructor for the System class. Initializes the system with the given Hamiltonian, wave function,
        // Monte Carlo solver, and list of particles.
        //
        // Input:   std::unique_ptr<class Hamiltonian> hamiltonian - Hamiltonian of the system;
        //          std::unique_ptr<class WaveFunction> waveFunction - wave function of the system;
        //          std::unique_ptr<class MonteCarlo> solver - solver for the Monte Carlo simulation;
        //          std::vector<std::unique_ptr<class Particle>> particles - list of particles.
System::System(
        std::unique_ptr<class Hamiltonian> hamiltonian,
        std::unique_ptr<class WaveFunction> waveFunction,
        std::unique_ptr<class MonteCarlo> solver,
        std::vector<std::unique_ptr<class Particle>> particles)
{
    m_numberOfParticles = particles.size();
    m_numberOfDimensions = particles[0]->getNumberOfDimensions();
    m_hamiltonian = std::move(hamiltonian);
    m_waveFunction = std::move(waveFunction);
    m_solver = std::move(solver);
    m_particles = std::move(particles);
}


// Performs equilibration steps for the system, adjusting the state based on the given step length
        // and number of equilibration steps.
        //
        // Input:   double stepLength - step size for each Monte Carlo step;
        //          unsigned int numberOfEquilibrationSteps - number of equilibration steps to perform.
        //
        // Output:  unsigned int - number of equilibration steps performed.
unsigned int System::runEquilibrationSteps(
        double stepLength,
        unsigned int numberOfEquilibrationSteps)
{
    unsigned int acceptedSteps = 0;

    // Run the specified number of burn-in VMC iterations.
    for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++)
    {
        acceptedSteps += m_solver->step(stepLength, *m_waveFunction, m_particles);
    }

    return acceptedSteps;
}


// Runs Metropolis steps for the system and returns a sampler object for collecting data during the simulation.
        //
        // Input:   double stepLength - step size for each Monte Carlo step;
        //          unsigned int numberOfMetropolisSteps - number of Metropolis steps to perform.
        //
        // Output:  std::unique_ptr<class Sampler> - a sampler for collecting data during the simulation.
std::unique_ptr<class Sampler> System::runMetropolisSteps(
        double stepLength,
        unsigned int numberOfMetropolisSteps)
{
    auto sampler = std::make_unique<Sampler>(
            m_numberOfParticles,
            m_numberOfDimensions,
            stepLength,
            numberOfMetropolisSteps);

    // Run the specified number of VMC iterations and sample observables.
    for (unsigned int i = 0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = m_solver->step(stepLength, *m_waveFunction, m_particles);

        sampler->sample(acceptedStep, this);
    }

    // Compute the averages of all sampled observables, including energy and other variables.
    sampler->computeAverages();

    return sampler;
}


// Computes the local energy of the system, which is required for variational Monte Carlo simulations.
//
// Output:  double - computed local energy of the system.
double System::computeLocalEnergy()
{
    return m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles);
}


// Computes the distance between particles i and j in the system.
//
// Input:   int i - index of the first particle;
//          int j - index of the second particle.
//
// Output:  double - distance between particles i and j.
double System::computerij(int i, int j)
{
    return m_waveFunction -> r_ij(m_particles, i, j);
}


// Computes the sum of the square of the distances (r^2) for all particles in the system.
//
// Output:  double - sum of r^2 for all particles in the system.
double System::computer2()
{
    double r2 = 0;
    for(int i = 0; i < m_numberOfParticles; i++)
    {
        r2 += m_waveFunction -> r_squared(m_particles, i);
    }

    return r2;
}


// Helper function that provides access to the current wave function parameters.
//
// Output:  const std::vector<double>& - the current wave function parameters.
const std::vector<double>& System::getWaveFunctionParameters()
{
    return m_waveFunction->getParameters();
}