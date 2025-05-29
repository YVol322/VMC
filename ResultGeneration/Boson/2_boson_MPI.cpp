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
#include <mpi.h>


using namespace std;
using namespace std::chrono;


int main(int argc, char** argv)
{    
    int seed = 2025;

    double omega = 1.0;
    double alpha = 0.5;
	double stepLength = 1;
    double energy;
    int size, my_rank;

    unsigned int numberOfDimensions = 1;
    unsigned int numberOfParticles;
    unsigned int numberOfMetropolisSteps = 1e4;
    unsigned int numberOfEquilibrationSteps = 1e3;

    std::vector<int> numberOfDims = {1, 2, 3};
    std::vector<bool> useNumericalWF = {0, 1};
    std::vector<bool> useMetropolisHastings = {0, 1};

    std::string path = "Results/Tables/Bosons/Program2";
    std::filesystem::create_directories(path);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    seed *= (my_rank + 1);
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
	                auto rng = std::make_unique<Random>(seed);
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
	                	    numberOfEquilibrationSteps / size);


                    auto sampler = system->runMetropolisSteps(
                            stepLength,
                            numberOfMetropolisSteps / size);
                    auto stop = high_resolution_clock::now();

                    auto duration = duration_cast<std::chrono::duration<double>>(stop - start);

                    energy = sampler->getEnergy();
                    double total_energy = 0.0;
                    MPI_Reduce( &energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                    if (my_rank == 0)
                    {
                        double mean_energy = total_energy / size;
                        sampler->setEnergy(mean_energy);
                        sampler->setTime(duration.count());
                        sampler -> printOutputToTerminal(*system);
                        outFile << mean_energy << "\t" << duration.count() << "\t" << numberOfParticles << endl;
                    }
                }
            }
        }
    }
    MPI_Finalize();

    return 0;
}