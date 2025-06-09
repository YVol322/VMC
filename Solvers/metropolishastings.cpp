#include "metropolishastings.h"     // Include "metropolishastings" header file with declarations.


MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool MetropolisHastings::step(
    double stepLength,
    class WaveFunction& waveFunction,
    std::vector<std::unique_ptr<class Particle>>& particles
)
{
    const double dt = stepLength;
    const double sqrtDt = std::sqrt(dt);
    const double D = 0.5;  // Constant D in atomic units (a.u.).

    int nParticles = particles.size();
    int nDims = particles[0]->getNumberOfDimensions();

    // Compute the wave function for the current (old) state.
    double WFold = waveFunction.evaluate(particles);

    // Initialize a vector to store quantum forces for all particles in the old state.
    std::vector<std::vector<double>> forcesOld(nParticles, std::vector<double>(nDims));

    // Loop over all particles and compute their quantum forces in the old state.
    for (int i = 0; i < nParticles; ++i)
    {
        forcesOld[i] = waveFunction.quantumForce(particles, i);
    }

    // Initialize the bool variable acceptedAny, which tracks if any step was accepted.
    bool acceptedAny = false;

    for (int i = 0; i < nParticles; i++)
    {
        // Save the old position of particle i.
        std::vector<double> Rold = particles[i]->getPosition();

        // Initialize a vector to store the new position of particle i.
        std::vector<double> Rnew(nDims);

        for (int d = 0; d < nDims; d++)
        {
            // Compute drift based on quantum force and timestep.
            double drift = D * forcesOld[i][d] * dt;

            // Generate a normally distributed random number.
            double gauss = m_rng->nextGaussian(0.0, 1.0);

            // Compute the new position of particle i in dimension d.
            Rnew[d] = Rold[d] + drift + gauss * sqrtDt;

            // Move the particle i to the new position.
            particles[i] -> adjustPosition(Rnew[d] - Rold[d], d);
        }

        // Compute the wave function for the new state after particle movement.
        double WFnew = waveFunction.evaluate(particles);

        // Compute the new quantum force for particle i.
        auto forcesNew = waveFunction.quantumForce(particles, i);

        // Compute Green's function ratio for the move.
        double greens = waveFunction.GreensFunctionRatio(Rnew, Rold, dt, forcesOld[i], forcesNew);

        // Compute the ratio of wave functions (squared).
        double psi2Ratio = (WFnew * WFnew) / (WFold * WFold);

        // Compute the acceptance probability.
        double acceptP = std::min(1.0, psi2Ratio * greens);

        // If a random number is less than or equal to the acceptance probability, accept the step.
        if (m_rng->nextDouble() <= acceptP)
        {
            // Update the old wave function value and quantum forces for particle i.
            WFold = WFnew;
            forcesOld[i] = std::move(forcesNew);

            // Mark that at least one step was accepted.
            acceptedAny = true;
        }
        else // Reject the step and move particle i back to its original position.
        {
            for (int d = 0; d < nDims; d++)
            {
                particles[i] -> adjustPosition(Rold[d] - Rnew[d], d);
            }
        }
    }

    // Return true if any step was accepted, otherwise false.
    return acceptedAny;
}