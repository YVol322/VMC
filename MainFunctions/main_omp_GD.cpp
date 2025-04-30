#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

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
#include <omp.h>

using namespace std;
using namespace std::chrono;

int main() {    
    int seed = 2025;

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 2;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e5;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e4;

    double omega = 1.0;
    double alpha = 0.5;

	int mode = 0;

	double stepLength = 1;
	double learning_rate = 1e-2;
	double stop_at = 1e-2;
	double max_iters = 1000;
	double iter = 0;
	double l2_norm = 4.1;

	int numberOfPairs = numberOfParticles * (numberOfParticles - 1) / 2;
    std::vector<double> beta(numberOfPairs, 0.42);
	std::vector<double> betaPJ(1, 0.446563);
    std::vector<double> grad_beta(numberOfPairs, 1);
    std::vector<double> grad_betaPJ(1, 1);

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
    double total_energyO = 0.0;
    double total_O = 0.0;
    std::vector<double> total_rij(numberOfPairs, 0.0);  // Initialize to 0
    std::vector<double> total_El_rij(numberOfPairs, 0.0);  // Initialize to 0

    while(iter < max_iters && l2_norm > stop_at) {
        // Reset local accumulations before starting each iteration
        std::fill(total_rij.begin(), total_rij.end(), 0.0);
        std::fill(total_El_rij.begin(), total_El_rij.end(), 0.0);
        total_energy = 0.0;
        total_O = 0.0;
        total_energyO = 0.0;

        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            int thread_seed = seed + thread_id;

            auto rng = std::make_unique<Random>(thread_seed);
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
                numberOfEquilibrationSteps/n_threads
            );

            auto sampler = system->runMetropolisSteps(
                stepLength,
                numberOfMetropolisSteps/n_threads
            );

            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::duration<double>>(stop - start);
            sampler->setTime(duration.count());

            double energy = sampler->getEnergy();
            std::vector<double> rij = sampler->getrij();
            std::vector<double> Elrij = sampler->getEnergyrij();
            double O = sampler->getO().at(0);
            double ElO = sampler->getEnergyO().at(0);

            #pragma omp atomic
            total_energy += energy;

            #pragma omp atomic
            total_O += O;

            #pragma omp atomic
            total_energyO += ElO;

            // Accumulate values for rij and El_rij
            for (int i = 0; i < numberOfPairs; i++) {
                #pragma omp atomic
                total_rij.at(i) += rij.at(i);

                #pragma omp atomic
                total_El_rij.at(i) += Elrij.at(i);
            }
            #pragma omp single
            {
                if(l2_norm < stop_at * 1.1) sampler -> printOutputToTerminal(*system);
            }
        }

        // Calculate the correct mean values across threads
        double mean_energy = total_energy / n_threads;
        double mean_energyO = total_energyO / n_threads;
        double mean_O = total_O / n_threads;

        std::vector<double> mean_rij(numberOfPairs, 0.0);  // Reset mean_rij
        std::vector<double> mean_El_rij(numberOfPairs, 0.0);  // Reset mean_El_rij

        for (int i = 0; i < numberOfPairs; i++) {
            mean_rij.at(i) = total_rij.at(i) / n_threads;
            mean_El_rij.at(i) = total_El_rij.at(i) / n_threads;
        }

        // Perform GD update outside parallel region for safety
        if(mode == 0)
        {
            for (int i = 0; i < numberOfPairs; i++) {
                grad_beta.at(i) = 2 * (mean_El_rij.at(i) - mean_energy * mean_rij.at(i));
                beta.at(i) -= learning_rate * grad_beta.at(i);
            }

            l2_norm = 0.0;
            for (double g : grad_beta) {
                l2_norm += g * g;
            }
            l2_norm = sqrt(l2_norm);

            std::cout << "Iteration " << iter
                      << ", beta = " << beta[0]
                      << ", energy = " << mean_energy
                      << ", grad = " << l2_norm << std::endl;

            iter++;
        }
        else
        {
            grad_betaPJ.at(0) = 2 * (mean_energyO - mean_energy * mean_O);
            betaPJ.at(0) -= learning_rate * grad_betaPJ.at(0);

            l2_norm = 0.0;
            for (double g : grad_betaPJ) {
                l2_norm += g * g;
            }
            l2_norm = sqrt(l2_norm);

            std::cout << "Iteration " << iter
                      << ", betaPJ = " << betaPJ[0]
                      << ", energy = " << mean_energy
                      << ", grad = " << l2_norm << std::endl;

            iter++;
        }
    }

    if (mode == 0) {
        for (double b : beta) std::cout << b << std::endl;
    } else {
        for (double x : betaPJ) std::cout << x << std::endl;
    }

    return 0;
}
