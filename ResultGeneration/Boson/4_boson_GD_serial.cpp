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

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e4;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e3;

    double omega = 1.0;
    double alpha = 0.5;

	int mode = 1;

	double stepLength = 1e-1;
	double learning_rate = 1e-1;
	double stop_at = 1e-4;
	double max_iters = 1000;
	double iter = 0;



	std::string path = "Results/Tables/Bosons/Program4";
    std::filesystem::create_directories(path);

	std::string filename = path + "/MH_N_Grad_El_2D.cvv";

	ofstream outFile(filename);
	outFile << "N\tGrad\tEnergy" << endl;
		
	
	for(int numberOfParticles = 1; numberOfParticles < 11; numberOfParticles++)
	{
		double El = 0;
		double l2_norm = 4.1;

		while(iter < max_iters && l2_norm > stop_at)
		{
			auto rng = std::make_unique<Random>(seed);
			auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);
			auto system = std::make_unique<System>(
			std::make_unique<HarmonicOscillator>(omega, 0),
			std::make_unique<Boson>(alpha, numberOfParticles),
			std::make_unique<MetropolisHastings>(std::move(rng)),
			std::move(particles));


			auto start = high_resolution_clock::now();
			auto acceptedEquilibrationSteps = system->runEquilibrationSteps(
					stepLength,
					numberOfEquilibrationSteps);

			auto sampler = system->runMetropolisSteps(
					stepLength,
					numberOfMetropolisSteps);
			auto stop = high_resolution_clock::now();

			auto duration = duration_cast<std::chrono::duration<double>>(stop - start);

			sampler -> setTime(duration.count());
			El = sampler -> getEnergy();
			double O1 = sampler -> getO1alpha();
			double O2 = sampler -> getO2alpha();

			double grad = 2 * (O2 - O1 * El);
			alpha -= learning_rate * grad;

			sampler -> printOutputToTerminal(*system);
			l2_norm = grad;
			std::cout << grad << std::endl;
		}
		outFile << numberOfParticles << "\t" << l2_norm << "\t" << El << endl;
	}

    return 0;
}