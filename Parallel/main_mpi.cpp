#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
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
    // Seed for the random number generator
    int seed = 2025;

    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 15;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e6;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e5;
    double omega = 1.0; // Oscillator frequency.
    double alpha = 0.5; // Variational parameter.
    double stepLength = 1; // Metropolis step length.

    int size, my_rank;
    double energy;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    seed *= (my_rank + 1);


    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);
    // Initialize particles
    auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);



    // Construct a unique pointer to a new System
    auto system = std::make_unique<System>(
            // Construct unique_ptr to Hamiltonian
            std::make_unique<HarmonicOscillator>(omega),
            // Construct unique_ptr to wave function
            std::make_unique<SimpleGaussian>(alpha),
            // Construct unique_ptr to solver, and move rng
            std::make_unique<Metropolis>(std::move(rng)),
            // Move the vector of particles to system
            std::move(particles));


    auto start = high_resolution_clock::now();


    // Run steps to equilibrate particles
    auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
            stepLength,
            numberOfEquilibrationSteps/size);

    // Run the Metropolis algorithm
    auto sampler = system->runMetropolisSteps(
            stepLength,
            numberOfMetropolisSteps/size);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    energy = sampler -> getEnergy();

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

        sampler -> setEnergy(mean_energy);
        sampler -> setTime(duration.count());
        sampler -> printOutputToTerminal(*system);

        delete[] all_energies;
        cout << " MPI" << endl;
    }

    MPI_Finalize();


    return 0;
}