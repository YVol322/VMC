#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
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
#include "Solvers/metropolishastings.h"
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

    double omega = 0.5;
    double alpha = 0.5;

	int mode = 1;
    bool useMetropolisHastings = 1;

	std::vector<double> stepLength_arr = {1};
	std::vector<double> learning_rate_arr = {1e-1};
	double stop_at = 1e-5;
	double max_iters = 100;
	double iter = 0;
	double l2_norm = 4.1;

	int numberOfPairs = numberOfParticles * (numberOfParticles - 1) / 2;
    std::vector<double> beta0(numberOfPairs, 0.381);
	std::vector<double> betaPJ0(1, 0.1);
    std::vector<double> grad_beta(numberOfPairs, 1);
    std::vector<double> grad_betaPJ(1, 1);

    std::string path = "Results/Tables/Fermions/Interacting/Program7";
    std::filesystem::create_directories(path);

    int n_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            n_threads = omp_get_num_threads();
            cout << "Number of threads: " << n_threads << endl;
        }
    }

    double total_El = 0.0;
    double total_O1Pade = 0.0;
    double total_O2Pade = 0.0;
    std::vector<double> total_O1Jastrow(numberOfPairs, 0.0);
    std::vector<double> total_O2Jastrow(numberOfPairs, 0.0);

    for(double stepLength: stepLength_arr)
    {
        for(double learning_rate: learning_rate_arr)
        {
            std::string filename = path
                + "/OMP_omega="
                + std::to_string(omega);

            if(useMetropolisHastings) filename += "_MH";
            if(mode) filename += "_PJ";

            std::ostringstream s;
            s << learning_rate;
            std::string eta_str = s.str(); 

            std::ostringstream ss;
            ss << stepLength;
            std::string h_str = ss.str(); 

            filename += "_eta=";
            filename += eta_str;

            filename += "_h=";
            filename += h_str;

            filename += ".cvv";


            ofstream outFile(filename);
            outFile << "Iter\tEnergy\tGrad\tBeta" << endl;


            l2_norm = 3;
            std::vector<double> beta = beta0;
            std::vector<double> betaPJ = betaPJ0;
            iter = 0;
            while(iter < max_iters && l2_norm > stop_at)
            {
                std::fill(total_O1Jastrow.begin(), total_O1Jastrow.end(), 0.0);
                std::fill(total_O2Jastrow.begin(), total_O2Jastrow.end(), 0.0);
                total_El = 0.0;
                total_O1Pade = 0.0;
                total_O2Pade = 0.0;

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
                        std::make_unique<HarmonicOscillator>(omega, 1),
                        std::make_unique<Fermion>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles, omega),
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
                    sampler->setTime(duration.count());

                    double energy = sampler -> getEnergy();
                    std::vector<double> O1Jastrow = sampler -> getO1Jastow();
                    std::vector<double> O2Jastrow = sampler -> getO2Jastow();
                    double O1Pade = sampler -> getO1Pade();
                    double O2Pade = sampler -> getO2Pade();

                    #pragma omp atomic
                    total_El += energy;

                    #pragma omp atomic
                    total_O1Pade += O1Pade;

                    #pragma omp atomic
                    total_O2Pade += O2Pade;

                    for (int i = 0; i < numberOfPairs; i++)
                    {
                        #pragma omp atomic
                        total_O1Jastrow[i] += O1Jastrow[i];

                        #pragma omp atomic
                        total_O2Jastrow[i] += O2Jastrow[i];
                    }
                }

                double mean_El = total_El / n_threads;
                double mean_O1Pade = total_O1Pade / n_threads;
                double mean_O2Pade = total_O2Pade / n_threads;

                std::vector<double> mean_O1Jastrow(numberOfPairs, 0.0);
                std::vector<double> mean_O2Jastrow(numberOfPairs, 0.0);

                for (int i = 0; i < numberOfPairs; i++)
                {
                    mean_O1Jastrow[i] = total_O1Jastrow[i]/ n_threads;
                    mean_O2Jastrow[i] = total_O2Jastrow[i]/ n_threads;
                }

                if(mode == 0)
                {
                    for (int i = 0; i < numberOfPairs; i++)
                    {
                        grad_beta[i] = 2 * (mean_O2Jastrow[i] - mean_El * mean_O1Jastrow[i]);
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
                            << ", energy = " << mean_El
                            << ", grad = " << l2_norm << std::endl;

                    iter++;
                }
                else
                {
                    grad_betaPJ[0] = 2 * (mean_O2Pade - mean_El * mean_O1Pade);
                    betaPJ[0] -= learning_rate * grad_betaPJ[0];

                    l2_norm = 0.0;
                    for (double g : grad_betaPJ)
                    {
                        l2_norm += g * g;
                    }
                    l2_norm = sqrt(l2_norm);

                    std::cout << "Iteration " << iter
                            << ", betaPJ = " << betaPJ[0]
                            << ", energy = " << mean_El
                            << ", grad = " << l2_norm << std::endl;

                    iter++;
                }
                outFile << iter << "\t" << total_El / n_threads << "\t" << l2_norm;
                if(mode == 0) outFile << "\t" << beta[0] << endl;
                else outFile << "\t" << betaPJ[0] << endl;
            }

            total_El = 0.0;
            unsigned int numberOfMetropolisSteps = (unsigned int) 1e7;
            unsigned int numberOfEquilibrationSteps = (unsigned int) 1e6;

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
                std::make_unique<HarmonicOscillator>(omega, 1),
                std::make_unique<Fermion>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles, omega),
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

                double energy = sampler -> getEnergy();

                #pragma omp atomic
                total_El += energy;
            }
            iter++;
            outFile << iter << "\t" << total_El / n_threads << "\t" << l2_norm;
            std::cout << total_El / n_threads << std::endl;
            
        }
    }

    return 0;
}