#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <mpi.h>
#include "system.h"
#include "WaveFunctions/fermionnumerical.h"
#include "WaveFunctions/fermion.h"
#include "WaveFunctions/boson.h"
#include "WaveFunctions/bosonnumerical.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Solvers/metropolishastings.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"
#include <autodiff/forward/dual.hpp>

using namespace std;
using namespace std::chrono;

int main(int argc, char** argv) {    
    int seed = 2025;

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 12;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e2;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e1;

    double omega = 1.0;
    double alpha = 0.5;

	int mode = 0;

	double stepLength = 1e-2;
	double learning_rate = 1e-3;
	double stop_at = 1e-1;
	double max_iters = 1000;
	double iter = 0;
	double l2_norm = 4.1;

	int numberOfPairs = numberOfParticles * (numberOfParticles - 1) / 2;
    std::vector<double> beta(numberOfPairs, 0.2);
	std::vector<double> betaPJ(1, 0.47);
    std::vector<double> grad_beta(numberOfPairs, 1);
    std::vector<double> grad_betaPJ(1, 1);

    double energy;
    double O1Pade;
    double O2Pade;
    std::vector<double> O1Jastrow(numberOfPairs, 0);
    std::vector<double> O2Jastrow(numberOfPairs, 0);

    int size, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    seed *= (my_rank + 1);

    while(iter < max_iters && l2_norm > stop_at)
    {
        auto rng = std::make_unique<Random>(seed);
        auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

        auto system = std::make_unique<System>(
            std::make_unique<HarmonicOscillator>(omega, 1),
            std::make_unique<Fermion>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles),
            std::make_unique<MetropolisHastings>(std::move(rng)),
            std::move(particles)
        );

        auto start = high_resolution_clock::now();
        auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
            stepLength,
            numberOfEquilibrationSteps / size
        );

        auto sampler = system->runMetropolisSteps(
            stepLength,
            numberOfMetropolisSteps / size
        );

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<std::chrono::duration<double>>(stop - start);

        energy = sampler->getEnergy();
        double total_energy = 0.0;
        MPI_Reduce( &energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (my_rank == 0)
        {
            sampler -> setEnergy(total_energy / size);
            sampler -> setTime(duration.count());
        }

        if (mode == 0)
        {
            O1Jastrow = sampler -> getO1Jastow();
            O2Jastrow = sampler -> getO2Jastow();

            std::vector<double> total_O1Jastrow(numberOfPairs);
            std::vector<double> total_O2Jastrow(numberOfPairs);

            MPI_Allreduce(O1Jastrow.data(), total_O1Jastrow.data(), 
                          numberOfPairs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            MPI_Allreduce(O2Jastrow.data(), total_O2Jastrow.data(), 
                          numberOfPairs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            

            if (my_rank == 0)
            {
                for (int i = 0; i < numberOfPairs; ++i)
                {
                    grad_beta[i] = 2 * (total_O2Jastrow[i]/size - energy * total_O1Jastrow[i] / size);
                    beta[i] -= learning_rate * grad_beta[i];
                }

                l2_norm = 0.0;
                for (double g : grad_beta)
                {
                    l2_norm += g * g;
                }
                l2_norm = sqrt(l2_norm);

                std::cout << "Iteration " << iter
                          << ", beta = " << beta[0]
                          << ", energy = " << energy
                          << ", grad = " << l2_norm << std::endl;

                iter++;
            }

            MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(l2_norm < stop_at && my_rank == 0) sampler -> printOutputToTerminal(*system);
        }
        else
        {
            O1Pade = sampler -> getO1Pade();
            O2Pade = sampler -> getO2Pade();

            double total_O1Pade = 0;
            double total_O2Pade = 0;

            MPI_Allreduce(&O1Pade, &total_O1Pade, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&O2Pade, &total_O2Pade, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (my_rank == 0)
            {
                grad_betaPJ[0] = 2 * (total_O2Pade / size - energy * total_O1Pade / size);
                betaPJ[0] -= learning_rate * grad_betaPJ[0];

                l2_norm = 0.0;
                for (double g : grad_betaPJ)
                {
                    l2_norm += g * g;
                }
                l2_norm = sqrt(l2_norm);

                std::cout << "Iteration " << iter
                          << ", betaPJ = " << betaPJ[0]
                          << ", energy = " << energy
                          << ", grad = " << l2_norm << std::endl;

                iter++;
            }

            MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if(l2_norm < stop_at && my_rank == 0) sampler -> printOutputToTerminal(*system);
    }

    MPI_Finalize();

    return 0;
}