#include "metropolis.h"
#include <iostream>


Metropolis::Metropolis(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool Metropolis::step(
            double stepLength,
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
        )
{
    int n_particles = particles.size();
    int n_dims = (particles.back()) -> getNumberOfDimensions();

    double WFold, WFnew;
    std::vector<std::vector<double>> rand_doubles(n_particles, std::vector<double>(n_dims));

    WFold = waveFunction.evaluate(particles);

    for(int i = 0; i < n_particles; i++)
    {
        for(int j = 0; j < n_dims; j++)
        {
            rand_doubles[i][j] = m_rng->nextDouble();

            (*(particles.at(i))).adjustPosition(stepLength*(rand_doubles[i][j] - 0.5), j);
        }
    }

    WFnew = waveFunction.evaluate(particles);

    double accept = std::min(1.0, (WFnew * WFnew)/(WFold * WFold));

    if(m_rng -> nextDouble() > accept)
    {
        for(int i = 0; i < n_particles; i++)
        {
            for(int j = 0; j < n_dims; j++)
            {
                (*(particles.at(i))).adjustPosition(-stepLength*(rand_doubles[i][j] - 0.5), j);
            }
        }
        return false;
    }
    else return true;

}