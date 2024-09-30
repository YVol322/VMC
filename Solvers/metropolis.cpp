#include <memory>
#include <vector>

#include "metropolis.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "Math/random.h"


Metropolis::Metropolis(std::unique_ptr<class Random> rng)
    : MonteCarlo(std::move(rng))
{
}


bool Metropolis::step(
        double stepLength,
        class WaveFunction& waveFunction,
        std::vector<std::unique_ptr<class Particle>>& particles)
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

    if(m_rng -> nextDouble() > (WFnew * WFnew)/(WFold * WFold))
    {
        //#pragma omp parallel for collapse(2)
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
