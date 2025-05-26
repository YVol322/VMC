#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
#include "WaveFunctions/fermionNI.h"
#include "WaveFunctions/fermionNInumerical.h"
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
#include <omp.h>

using namespace std;
using namespace std::chrono;

int main()
{    
    int seed = 2025;

    double omega = 1.0;
    double alpha = 0.52;

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 12;
    unsigned int numberOfMetropolisSteps = 1e3;
    unsigned int numberOfEquilibrationSteps = 1e2;

    bool useMetropolisHastings = 0;
    bool useNum = 0;

	double stepLength = 1e-1;
	double learning_rate = 1e-2;
	double stop_at = 1e-3;
	double l2_norm = 4.1;


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
    double total_O1 = 0.0;
    double total_O2 = 0.0;

    int max_iters = 10000;
    int iter = 0;

    double elapsed_time;

    while(iter < max_iters && l2_norm > stop_at)
    {
        total_energy = 0.0;
        total_O1 = 0.0;
        total_O2 = 0.0;

        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            int thread_seed = seed + thread_id;

            auto rng = std::make_unique<Random>(thread_seed);
            auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

            unique_ptr<MonteCarlo> solver;
            unique_ptr<WaveFunction> WF;

            if(useMetropolisHastings) solver = std::make_unique<MetropolisHastings>(std::move(rng));
            else solver = std::make_unique<Metropolis>(std::move(rng));

            if(useNum) WF = std::make_unique<FermionNInumerical>(alpha, numberOfParticles);
            else WF = std::make_unique<FermionNI>(alpha, numberOfParticles);

            auto system = std::make_unique<System>(
                std::make_unique<HarmonicOscillator>(omega),
                std::move(WF),
                std::move(solver),
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
            elapsed_time = duration.count();
            sampler->setTime(elapsed_time);

            double energy = sampler -> getEnergy();
            double O1 = sampler -> getO1alpha();
            double O2 = sampler -> getO2alpha();

            #pragma omp atomic
            total_energy += energy;

            #pragma omp atomic
            total_O1 += O1;

            #pragma omp atomic
            total_O2 += O2;

            #pragma omp critical
            elapsed_time = duration.count();
        }

        double mean_energy = total_energy / n_threads;
        double mean_O1 = total_O1 / n_threads;
        double mean_O2 = total_O2 / n_threads;

        double grad = 2 * (mean_O2 - mean_O1 * mean_energy);
        std::cout << "iter = " << iter << std::endl;
	    std::cout << "local en = "<< mean_energy << std::endl;
	    std::cout << "grad = "<< grad << std::endl;
	    std::cout << "alpha = "<< alpha << std::endl;
	    //std::cout << "time = "<< elapsed_time << std::endl;
	    std::cout << std::endl;
	    alpha += learning_rate * grad;

        l2_norm = abs(grad);

        iter++;
    }

    return 0;
}
