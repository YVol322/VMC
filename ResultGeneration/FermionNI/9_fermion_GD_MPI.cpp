#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <mpi.h>
#include "system.h"
#include "WaveFunctions/fermionNInumerical.h"
#include "WaveFunctions/fermionNI.h"
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
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e3;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e2;

    double omega = 1.0;
    double alpha = 0.6;

	int mode = 0;

	double stepLength = 5e-2;
	double learning_rate = 1e-2;
	double stop_at = 1e-3;
	double max_iters = 1000;
	double iter = 0;
	double l2_norm = 4.1;

    double energy = 0;
    double O1 = 0;
    double O2 = 0;
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
            std::make_unique<HarmonicOscillator>(omega, 0),
            std::make_unique<FermionNInumerical>(alpha, numberOfParticles),
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

        O1 = sampler->getO1alpha();
        double total_O1 = 0.0;
        MPI_Reduce( &O1, &total_O1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        O2 = sampler->getO2alpha();
        double total_O2 = 0.0;
        MPI_Reduce( &O2, &total_O2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0)
        {
            double mean_energy = total_energy / size;
            double mean_O1 = total_O1 / size;
            double mean_O2 = total_O2 / size;
            sampler -> setEnergy(mean_energy);
            sampler -> setTime(duration.count());

            double grad = 2 * (mean_O2 - mean_O1 * mean_energy);

            alpha += grad * learning_rate;

            l2_norm = abs(grad);

            std::cout << "iter = " << iter << std::endl;
            std::cout << "local en = "<< total_energy / size << std::endl;
            std::cout << "grad = "<< grad << std::endl;
            std::cout << "alpha = "<< alpha << std::endl;
            std::cout << "time = "<< duration.count() << std::endl;
            iter++;
        }
        MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}