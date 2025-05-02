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
    unsigned int numberOfParticles = 6;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e3;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e2;

    double omega = 1.0;
    double alpha = 0.5;

	int mode = 1;

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

    int size, my_rank;
    double energy;
    double O;
    double energyO;
    std::vector<double> rij(numberOfPairs, 0);
    std::vector<double> Erij(numberOfPairs, 0);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    seed *= (my_rank + 1);

    while(iter < max_iters && l2_norm > stop_at)
    {
        // Initialize random engine for each process
        auto rng = std::make_unique<Random>(seed);
        auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

        auto system = std::make_unique<System>(
            std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<Fermion>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles),
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
        if (my_rank == 0) {
            // Root process allocates space to store all energies
            all_energies = new double[size];
        }

        MPI_Gather(&energy, 1, MPI_DOUBLE, all_energies, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (my_rank == 0) {
            double sum_energy = 0.0;
            for (int i = 0; i < size; ++i) {
                sum_energy += all_energies[i];
            }
            double mean_energy = sum_energy / size;
            sampler->setEnergy(mean_energy);
            sampler->setTime(duration.count());
        }

        // Mode 0 for beta gradients
        if (mode == 0)
        {
            rij = sampler->getrij();
            Erij = sampler->getEnergyrij();

            // Allocate memory for all processes to send their data
            double *all_rijs = nullptr;
            double *all_energyrijs = nullptr;
            if (my_rank == 0) {
                // Root process allocates space to store all rijs and energyrijs
                all_rijs = new double[size * rij.size()];
                all_energyrijs = new double[size * Erij.size()];
            }

            // Gather rij and Erij from all processes to the root process
            MPI_Gather(rij.data(), rij.size(), MPI_DOUBLE, all_rijs, rij.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(Erij.data(), Erij.size(), MPI_DOUBLE, all_energyrijs, Erij.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (my_rank == 0) {
                double sum_rij = 0.0;
                double sum_El_rij = 0.0;
                // Compute averages for rijs and energyrijs
                for (int i = 0; i < size; ++i) {
                    sum_rij += all_rijs[i];  // Add all rijs for this pair across all processes
                    sum_El_rij += all_energyrijs[i];  // Add all Erij for this pair across all processes
                }

                // Compute means for all rijs and energyrijs
                double mean_rij = sum_rij / size;
                double mean_El_rij = sum_El_rij / size;

                // Update gradients
                for (int i = 0; i < numberOfPairs; ++i) {
                    grad_beta[i] = 2 * (mean_El_rij - energy * mean_rij);
                    beta[i] -= learning_rate * grad_beta[i];
                }

                // Compute L2 norm
                l2_norm = 0.0;
                for (double g : grad_beta) {
                    l2_norm += g * g;
                }
                l2_norm = sqrt(l2_norm);

                std::cout << "Iteration " << iter
                          << ", beta = " << beta[0]
                          << ", energy = " << energy
                          << ", grad = " << l2_norm << std::endl;

                iter++;

                delete[] all_rijs;
                delete[] all_energyrijs;
            }

            MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(l2_norm < stop_at && my_rank == 0) sampler -> printOutputToTerminal(*system);
        }
        else
        {
            // Mode 1 for updating betaPJ
            O = sampler->getO().at(0);
            energyO = sampler->getEnergyO().at(0);

            // Allocate memory for all processes to send their data
            double *all_Os = nullptr;
            double *all_energyOs = nullptr;
            if (my_rank == 0)
            {
                // Root process allocates space to store all Os and energyOs
                all_Os = new double[size];
                all_energyOs = new double[size];
            }

            // Gather O and energyO from all processes to the root process
            MPI_Gather(&O, 1, MPI_DOUBLE, all_Os, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(&energyO, 1, MPI_DOUBLE, all_energyOs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (my_rank == 0)
            {
                double sum_O = 0.0;
                double sum_energyO = 0.0;
                // Compute averages for O and energyO
                for (int i = 0; i < size; ++i)
                {
                    sum_O += all_Os[i];  // Add all O values across all processes
                    sum_energyO += all_energyOs[i];  // Add all energyO values across all processes
                }

                // Compute means for O and energyO
                double mean_O = sum_O / size;
                double mean_energyO = sum_energyO / size;

                // Update grad_betaPJ
                grad_betaPJ[0] = 2 * (mean_energyO - energy * mean_O);
                betaPJ[0] -= learning_rate * grad_betaPJ[0];

                // Compute L2 norm for grad_betaPJ
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

                delete[] all_Os;
                delete[] all_energyOs;
            }

            MPI_Bcast(&l2_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if(l2_norm < stop_at && my_rank == 0) sampler -> printOutputToTerminal(*system);
    }

    MPI_Finalize();

    return 0;
}