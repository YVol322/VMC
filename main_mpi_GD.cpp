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
    double alpha = 0.5;

	int mode = 0;

	double stepLength = 0.5;
	double learning_rate = 1e-5;
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
            std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<FermionNumerical>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles),
            std::make_unique<Metropolis>(std::move(rng)),
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
        double *all_energies = nullptr;
        if (my_rank == 0)
        {
            all_energies = new double[size];
        }

        MPI_Gather(&energy, 1, MPI_DOUBLE, all_energies, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (my_rank == 0)
        {
            double sum_energy = 0.0;

            for (int i = 0; i < size; ++i)
            {
                sum_energy += all_energies[i];
            }

            double mean_energy = sum_energy / size;
            sampler -> setEnergy(mean_energy);
            sampler -> setTime(duration.count());
        }

        if (mode == 0)
        {
            O1Jastrow = sampler -> getO1Jastow();
            O2Jastrow = sampler -> getO2Jastow();


            double *all_O1Js = nullptr;
            double *all_O2Js = nullptr;
            if (my_rank == 0)
            {
                all_O1Js = new double[size * O1Jastrow.size()];
                all_O2Js = new double[size * O2Jastrow.size()];
            }

            MPI_Gather(O1Jastrow.data(), O1Jastrow.size(), MPI_DOUBLE, all_O1Js, O1Jastrow.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(O2Jastrow.data(), O2Jastrow.size(), MPI_DOUBLE, all_O2Js, O2Jastrow.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (my_rank == 0)
            {
                double sum_O1J = 0.0;
                double sum_O2J = 0.0;

                for (int i = 0; i < size; ++i)
                {
                    sum_O1J += all_O1Js[i];
                    sum_O2J += all_O2Js[i];
                }

                double mean_O1 = sum_O1J / size;
                double mean_O2 = sum_O2J / size;

                for (int i = 0; i < numberOfPairs; ++i)
                {
                    grad_beta[i] = 2 * (mean_O2 - energy * mean_O1);
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

                delete[] all_O1Js;
                delete[] all_O2Js;
            }

            MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(l2_norm < stop_at && my_rank == 0) sampler -> printOutputToTerminal(*system);
        }
        else
        {
            O1Pade = sampler -> getO1Pade();
            O2Pade = sampler -> getO2Pade();

            double *all_O1Ps = nullptr;
            double *all_O2Ps = nullptr;
            if (my_rank == 0)
            {
                all_O1Ps = new double[size];
                all_O2Ps = new double[size];
            }

            MPI_Gather(&O1Pade, 1, MPI_DOUBLE, all_O1Ps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(&O2Pade, 1, MPI_DOUBLE, all_O2Ps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (my_rank == 0)
            {
                double sum_O1P = 0.0;
                double sum_O2P = 0.0;
                for (int i = 0; i < size; ++i)
                {
                    sum_O1P += all_O1Ps[i];
                    sum_O2P += all_O2Ps[i];
                }

                double mean_O1P = sum_O1P / size;
                double mean_O2P = sum_O2P / size;

                grad_betaPJ[0] = 2 * (mean_O2P - energy * mean_O1P);
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

                delete[] all_O1Ps;
                delete[] all_O2Ps;
            }

            MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if(l2_norm < stop_at && my_rank == 0) sampler -> printOutputToTerminal(*system);
    }

    MPI_Finalize();

    return 0;
}