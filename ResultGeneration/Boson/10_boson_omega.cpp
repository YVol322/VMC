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
    unsigned int numberOfDimensions = 3;
    unsigned int numberOfMetropolisSteps = 1e4;
    unsigned int numberOfEquilibrationSteps = 1e3;

    std::vector<double> omegas = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<bool> useNumericalWF = {0};
    std::vector<bool> useMetropolisHastings = {0};

    std::string path = "Results/Tables/Bosons/Program7";
    std::filesystem::create_directories(path);
    std::string filename = path
                    + "/"
                    + std::to_string(numberOfDimensions)
                    + "D_";

                filename += "Energy_omega";

                filename += std::to_string(omega);
                filename += ".cvv";

                ofstream outFile(filename);
                outFile << "Energy\tOmega" << endl;

    for(bool MH: useMetropolisHastings)
    {
        for(bool NUM: useNumericalWF)
        {
            for(double omega: omegas)
            {

                for (int i = 3; i < 4; i++)
                {
                    numberOfParticles = i;

                    auto rng = make_unique<Random>(seed);
                    auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);

                    unique_ptr<WaveFunction> wf;
                    unique_ptr<MonteCarlo> solver;

                    if(NUM) wf = make_unique<BosonNumerical>(alpha, numberOfParticles, omega);
                    else wf = make_unique<Boson>(alpha, numberOfParticles, omega);

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

                    outFile << E_L << "\t" << omega << endl;
                }
            }
            outFile.close();
        }
    }
    return 0;
}