#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/slaterdeterminant.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"

#include <omp.h>

using namespace std;
using namespace std::chrono;


int main() {
    // Seed for the random number generator
    int seed = 2024;

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 6;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e5;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e4;
    double omega = 1.0; // Oscillator frequency.
    double alpha = 0.5; // Variational parameter.
    double stepLength = 1; // Metropolis step length.

    int n_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            n_threads = omp_get_num_threads();
            cout << "Number of threads: " << n_threads << endl;
        }
    }

    double total_energy = 0.0;

    auto start = high_resolution_clock::now();

    // Parallel region
    #pragma omp parallel reduction(+:total_energy)
    {
        int thread_id = omp_get_thread_num();
        // Create a unique seed for each thread based on the base seed
        int thread_seed = seed + thread_id;

        // Create a unique RNG for each thread
        auto rng = std::make_unique<Random>(thread_seed);

        // Initialize particles for this thread's system
        auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

        // Create a unique system for this thread
        auto system = std::make_unique<System>(
            std::make_unique<HarmonicOscillator>(omega),
            std::make_unique<SlaterDeterminant>(alpha),
            std::make_unique<Metropolis>(std::move(rng)),
            std::move(particles)
        );

        // Run equilibration steps
        auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
            stepLength,
            numberOfEquilibrationSteps/n_threads
        );

        // Run Metropolis steps
        auto sampler = system->runMetropolisSteps(
            stepLength,
            numberOfMetropolisSteps/n_threads
        );

        // (Optional) Store or process results from this thread's system
        // For example, you can use a reduction or critical section to aggregate results

        double energy = sampler->getEnergy();
        total_energy += energy;
        #pragma omp single
        {
            sampler->printOutputToTerminal(*system);
        }
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    double mean_energy = total_energy / n_threads;

    cout << "Parallel execution time: " << duration.count() << " seconds" << endl;
    cout << "Mean energy across all systems: " << mean_energy << endl;

    // Output information from the simulation
    //sampler->printOutputToTerminal(*system);

    return 0;
}