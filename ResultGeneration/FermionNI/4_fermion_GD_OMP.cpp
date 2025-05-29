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
    double alpha = 0.5;

    unsigned int numberOfDimensions = 2;
    std::vector<int> numberOfParticles_arr = {2, 6, 12};
    unsigned int numberOfMetropolisSteps = 1e3;
    unsigned int numberOfEquilibrationSteps = 1e2;

    std::vector<bool> useMetropolisHastings_arr = {0,1};

	double stepLength = 1e-1;
	double learning_rate = 1e-2;
	double stop_at = 1e-3;
	double l2_norm = 4.1;
    int max_iters = 3000;

    double grad = 0;
    double elapsed_time;

    int n_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            n_threads = omp_get_num_threads();
            cout << "Number of threads: " << n_threads << endl;
        }
    }

    std::string path = "Results/Tables/Fermions/NonInteracting/Program4";
    std::filesystem::create_directories(path);

    std::string filename = path
        + "/OMP_GD.cvv";

    ofstream outFile(filename);
    outFile << "Energy\tTime\tGrad" << endl;

    for(int numberOfParticles: numberOfParticles_arr)
    {
        for(bool useMetropolisHastings: useMetropolisHastings_arr)
        {
            double total_energy = 0.0;
            double total_O1 = 0.0;
            double total_O2 = 0.0;

            int iter = 0;
            l2_norm = 2;
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

                    if(useMetropolisHastings) solver = std::make_unique<MetropolisHastings>(std::move(rng));
                    else solver = std::make_unique<Metropolis>(std::move(rng));

                    auto system = std::make_unique<System>(
                        std::make_unique<HarmonicOscillator>(omega, 0),
                        std::make_unique<FermionNI>(alpha, numberOfParticles),
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

                grad = 2 * (mean_O2 - mean_O1 * mean_energy);
                std::cout << "iter = " << iter << std::endl;
                std::cout << "local en = "<< mean_energy << std::endl;
                std::cout << "grad = "<< grad << std::endl;
                std::cout << "alpha = "<< alpha << std::endl;
                std::cout << "time = "<< elapsed_time << std::endl;
                std::cout << std::endl;
                alpha += learning_rate * grad;

                l2_norm = abs(grad);

                iter++;
            }
            outFile << total_energy/n_threads << "\t" << elapsed_time << "\t" << grad << endl;
        }
    }

    return 0;
}
