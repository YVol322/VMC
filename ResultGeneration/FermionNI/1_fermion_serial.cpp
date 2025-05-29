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

using namespace std;
using namespace std::chrono;

int main()
{
    int seed = 2025;
    double omega = 1;
    double alpha = 0.5;
    double stepLength = 0.1;

    unsigned int numberOfParticles = 12;
    unsigned int numberOfDimensions = 2;
    unsigned int numberOfMetropolisSteps = 1e4;
    unsigned int numberOfEquilibrationSteps = 1e3;

    std::vector<bool> useNumericalWF = {0, 1};
    std::vector<bool> useMetropolisHastings = {0, 1};

    std::string path = "Results/Tables/Fermions/NonInteracting/Program1";
    std::filesystem::create_directories(path);

    std::string filename = path
                + "/N="
                + std::to_string(numberOfParticles)
                + ".cvv";

    ofstream outFile(filename);
    outFile << "Energy\tTime\tAccratio" << endl;

    for(bool MH: useMetropolisHastings)
    {
        for(bool NUM: useNumericalWF)
        {
            auto rng = make_unique<Random>(seed);
            auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

            unique_ptr<WaveFunction> wf;
            unique_ptr<MonteCarlo> solver;

            if(NUM) wf = make_unique<FermionNInumerical>(alpha, numberOfParticles);
            else wf = make_unique<FermionNI>(alpha, numberOfParticles);

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

            outFile << E_L << "\t" << timeElapsed << "\t" << accratio  << endl;

        }
    }
    outFile.close();

    return 0;
}