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
    double alpha = 0.5;
    double stepLength = 0.1;

    unsigned int numberOfParticles = 12;
    unsigned int numberOfDimensions = 2;
    unsigned int numberOfMetropolisSteps = 1e4;
    unsigned int numberOfEquilibrationSteps = 1e3;

    std::vector<double> omegas = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    std::vector<bool> useNumericalWF = {0};
    std::vector<bool> useMetropolisHastings = {0};

    std::string path = "Results/Tables/Fermions/NonInteracting/Program7";
    std::filesystem::create_directories(path);

    std::string filename = path
                + "/N="
                + std::to_string(numberOfParticles)
                + ".cvv";

    ofstream outFile(filename);
    outFile << "Energy\tOmega" << endl;

    for(double omega: omegas)
    {
        for(bool NUM: useNumericalWF)
        {
            auto rng = make_unique<Random>(seed);
            auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

            unique_ptr<WaveFunction> wf;
            unique_ptr<MonteCarlo> solver;

            if(NUM) wf = make_unique<FermionNInumerical>(alpha, numberOfParticles, omega);
            else wf = make_unique<FermionNI>(alpha, numberOfParticles, omega);
            
            solver = make_unique<Metropolis>(move(rng));

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

            outFile << E_L << "\t" << omega << endl;

        }
    }
    outFile.close();

    return 0;
}