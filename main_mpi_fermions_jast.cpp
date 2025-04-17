#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>

#include "system.h"
#include "WaveFunctions/fermionsjastrow.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

#include <mpi.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char** argv)
{
    int seed = 2025;
    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 6;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e4;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e3;
    double omega = 1.0;
    double alpha = 0.5;
    double stepLength = 1;

    int size, my_rank;
    double energy;

    int max_iterations = 1000;
    double learning_rate = 5e-4;
    double momentum = 0.9;
    double convergence_threshold = 0.9;
    double gradient_norm = 1.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    seed *= (my_rank + 1);

    std::vector<double> beta(15);
    std::vector<double> energyrij;
    std::vector<double> rij;
    auto rng2 = std::make_unique<Random>(seed);

    // Use good starting values
    if (my_rank == 0) {
        for (int i = 0; i < 15; ++i) {
            beta[i] = rng2 -> nextDouble() * 0.3;
        }
    }

    int iteration = 0;
    while (iteration < max_iterations && gradient_norm > convergence_threshold) {
        auto rng = std::make_unique<Random>(seed + iteration);
        auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

        auto system = std::make_unique<System>(
            std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<FermionsJastrow>(alpha, beta),
            std::make_unique<Metropolis>(std::move(rng)),
            std::move(particles));

        auto start = high_resolution_clock::now();

        system->runEquilibrationSteps(stepLength, numberOfEquilibrationSteps / size);
        auto sampler = system->runMetropolisSteps(stepLength, numberOfMetropolisSteps / size);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);

        energy = sampler->getEnergy();

        double* all_energies = nullptr;
        if (my_rank == 0) {
            all_energies = new double[size];
        }

        MPI_Gather(&energy, 1, MPI_DOUBLE, all_energies, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double mean_energy = energy;
        if (my_rank == 0) {
            double sum_energy = 0.0;
            for (int i = 0; i < size; ++i) {
                sum_energy += all_energies[i];
            }
            mean_energy = sum_energy / size;
            delete[] all_energies;
        }

        MPI_Bcast(&mean_energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        rij = sampler->getrij();
        energyrij = sampler->getEnergyrij();

        std::vector<double> gradient(15);
        std::vector<double> velocity(15);
        gradient_norm = 0.0;
        for (int k = 0; k < 15; ++k) {
            gradient[k] = 2.0 * (energyrij[k] - mean_energy * rij[k]);
            gradient_norm += gradient[k] * gradient[k];
        }
        gradient_norm = sqrt(gradient_norm);

        if (my_rank == 0) {
            for (int k = 0; k < 15; ++k) {
                velocity[k] = momentum * velocity[k] - learning_rate * gradient[k];
                beta[k] += velocity[k];
            }

            sampler->setEnergy(mean_energy);
            sampler->setTime(duration.count());

            cout << "SGD Iteration: " << iteration << endl;
            cout << "Energy: " << mean_energy << endl;
            cout << "Gradient norm: " << gradient_norm << endl;
            for (int k = 0; k < 15; ++k) {
                cout << "  beta[" << k << "] = " << beta[k] << " (grad = " << gradient[k] << ")" << endl;
            }
        }

        iteration++;
    }

    if (my_rank == 0) {
        if (gradient_norm <= convergence_threshold) {
            cout << "Converged after " << iteration << " iterations." << endl;
        } else {
            cout << "Max iterations reached without convergence." << endl;
        }
    }

    MPI_Finalize();
    return 0;
}