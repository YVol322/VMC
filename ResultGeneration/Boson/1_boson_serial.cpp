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

using namespace std;
using namespace std::chrono;

int main()
{
    int seed = 2025;
    double omega = 1;
    double alpha = 0.5;
    double stepLength = 1;

    unsigned int numberOfParticles;
    unsigned int numberOfMetropolisSteps = 1e4;
    unsigned int numberOfEquilibrationSteps = 1e3;

    std::vector<int> numberOfDims = {1, 2, 3};
    std::vector<bool> useNumericalWF = {0, 1};
    std::vector<bool> useMetropolisHastings = {0, 1};

    std::string path = "Results/Tables/Bosons/Program1";
    std::filesystem::create_directories(path);

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

                for (int i = 1; i < 11; i++)
                {
                    numberOfParticles = i;

                    auto rng = make_unique<Random>(seed);
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
                        numberOfEquilibrationSteps
                    );

                    auto sampler = system->runMetropolisSteps(
                        stepLength,
                        numberOfMetropolisSteps
                    );

                    auto stop = high_resolution_clock::now();
                    duration<double> elapsed = stop - start;
                    double timeElapsed = elapsed.count();
                    double accratio = double(acceptedEquilibrationSteps) / double(numberOfEquilibrationSteps);
                    sampler -> setTime(timeElapsed);


                    sampler -> printOutputToTerminal(*system);

                    double E_L = sampler->getEnergy();

                    outFile << E_L << "\t" << timeElapsed << "\t" << numberOfParticles << "\t" << accratio << endl;
                }


                outFile.close();
            }
        }
    }
    return 0;
}