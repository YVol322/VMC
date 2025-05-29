#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
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
    double alpha = 0.5;
	double stepLength = 1;

    unsigned int numberOfParticles;
    unsigned int numberOfMetropolisSteps = 1e4;
    unsigned int numberOfEquilibrationSteps = 1e3;

    std::vector<int> numberOfDims = {1, 2, 3};
    std::vector<bool> useNumericalWF = {0, 1};
    std::vector<bool> useMetropolisHastings = {0, 1};

    std::string path = "Results/Tables/Bosons/Program3";
    std::filesystem::create_directories(path);

    double elapsed_time, total_energy = 0.0;

    int n_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            n_threads = omp_get_num_threads();
            cout << "Number of threads: " << n_threads << endl;
        }
    }

    for(bool MH: useMetropolisHastings)
    {
        for(bool NUM: useNumericalWF)
        {
            for(int numberOfDimensions: numberOfDims)
            {
                std::string filename = path
                    + "/"
                    + std::to_string(numberOfDimensions)
                    + "D_";

                if(NUM) filename += "Num_";
                if(MH) filename += "MH_";
                filename += "Energy_Time_AccRate.cvv";

                ofstream outFile(filename);
                outFile << "Energy\tTime\tN\tAccRatio" << endl;
                for(int i = 1; i < 11; i++)
                {
                    numberOfParticles = i;
                    total_energy = 0;
                    #pragma omp parallel
                    {
                        int thread_id = omp_get_thread_num();
                        int thread_seed = seed + thread_id;

                        auto rng = std::make_unique<Random>(thread_seed);
                        auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

                        unique_ptr<WaveFunction> wf;
                        unique_ptr<MonteCarlo> solver;
                        if(NUM) wf = make_unique<BosonNumerical>(alpha, numberOfParticles);
                        else wf = make_unique<Boson>(alpha, numberOfParticles);

                        if(MH) solver = make_unique<MetropolisHastings>(move(rng));
                        else solver = make_unique<Metropolis>(move(rng));

                        auto system = make_unique<System>(
                            make_unique<HarmonicOscillator>(omega, 0),
                            move(wf),
                            move(solver),
                            move(particles)
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

                        double energy = sampler -> getEnergy();

                        #pragma omp atomic
                        total_energy += energy;

                        #pragma omp barrier

                        #pragma omp master
                        {
                            elapsed_time = duration.count();
                            sampler -> setEnergy(total_energy/n_threads);
                            sampler -> setTime(elapsed_time);
                            sampler ->printOutputToTerminal(*system);
                        }
                    }
                    outFile << total_energy/n_threads << "\t" << elapsed_time << "\t" << numberOfParticles << endl;
                }
            }
        }
    }

return 0;
}