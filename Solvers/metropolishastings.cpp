#include "metropolishastings.h"
#include <iostream>


MetropolisHastings::MetropolisHastings(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool MetropolisHastings::step(
    double stepLength,
    class WaveFunction& waveFunction,
    std::vector<std::unique_ptr<class Particle>>& particles
) {
    const double dt = stepLength;
    const double sqrtDt = std::sqrt(dt);
    const double D = 0.5;

    int nParticles = particles.size();
    int nDims = particles[0]->getNumberOfDimensions();

    double WFold = waveFunction.evaluate(particles);
    std::vector<std::vector<double>> forcesOld(nParticles, std::vector<double>(nDims));
    for (int i = 0; i < nParticles; ++i)
    {
        forcesOld[i] = waveFunction.quantumForce(particles, i);
    }

    bool acceptedAny = false;

    for (int i = 0; i < nParticles; i++)
    {
        std::vector<double> Rold = particles[i]->getPosition();

        std::vector<double> Rnew(nDims);
        for (int d = 0; d < nDims; d++)
        {
            double drift = D * forcesOld[i][d] * dt;
            double gauss = m_rng->nextGaussian(0.0, 1.0);
            Rnew[d] = Rold[d] + drift + gauss * sqrtDt;

            particles[i] -> adjustPosition(Rnew[d] - Rold[d], d);
        }

        double WFnew = waveFunction.evaluate(particles);
        auto   forcesNew = waveFunction.quantumForce(particles, i);

        double greens    = waveFunction.GreensFunctionRatio(Rnew, Rold, dt, forcesOld[i], forcesNew);
        double psi2Ratio = (WFnew * WFnew) / (WFold * WFold);
        double acceptP   = std::min(1.0, psi2Ratio * greens);

        if (m_rng->nextDouble() <= acceptP)
        {
            WFold = WFnew;
            forcesOld[i] = std::move(forcesNew);
            acceptedAny = true;
        }
        else
        {
            for (int d = 0; d < nDims; d++)
            {
                particles[i] -> adjustPosition(Rold[d] - Rnew[d], d);
            }
        }
    }

    return acceptedAny;
}
