#include <iostream>
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


using namespace std;
using namespace std::chrono;


int main()
{    
    int seed = 2025;

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 12;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e2;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e1;

    double omega = 1.0;
    double alpha = 0.5;

	int mode = 1;

	double stepLength = 1e-3;
	double learning_rate = 1e-3;
	double stop_at = 1e-2;
	double max_iters = 1000;
	double iter = 0;
	double l2_norm = 4.1;

	int numberOfPairs = numberOfParticles * (numberOfParticles - 1) / 2;
    std::vector<double> beta(numberOfPairs, 0.3);
	std::vector<double> betaPJ(1, 0.446563);
    std::vector<double> grad_beta(numberOfPairs, 1);
    std::vector<double> grad_betaPJ(1, 1);
	
	while(iter < max_iters && l2_norm > stop_at)
	{
		auto rng = std::make_unique<Random>(seed);
		auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);
		auto system = std::make_unique<System>(
        std::make_unique<HarmonicOscillator>(omega, 1),
        std::make_unique<FermionNumerical>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles),
		//std::make_unique<BosonNumerical>(alpha, betaPJ, mode, numberOfParticles),
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
		if(mode == 0)
		{
			double El = sampler -> getEnergy();
			std::vector<double> O1 = sampler -> getO1Jastow();
			std::vector<double> O2 = sampler -> getO2Jastow();
			for (int i = 0; i < numberOfPairs; i++)
			{
				grad_beta[i]= 2 * (O2[i] - O1[i] * El);

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
        	  << ", energy = " << El
        	  << ", grad = " << l2_norm << std::endl;

			iter++;

			if(l2_norm < stop_at)
			{
				sampler -> printOutputToTerminal(*system);
				return 0;
			}
		}
		else
		{
			double El = sampler -> getEnergy();
			double O1 = sampler -> getO1Pade();
			double O2 = sampler -> getO2Pade();
	
			grad_betaPJ[0] = 2 * (O2 - El * O1);
	
			betaPJ[0] -= learning_rate * grad_betaPJ[0];
	
			l2_norm = 0.0;
			for (double g : grad_betaPJ)
			{
			    l2_norm += g * g;
			}
			l2_norm = sqrt(l2_norm);

			std::cout << "Iteration " << iter
        	  << ", beta = " << betaPJ[0]
        	  << ", energy = " << El
        	  << ", grad = " << l2_norm << std::endl;

			iter++;
		}

		if(l2_norm < stop_at) sampler -> printOutputToTerminal(*system);
	}

	//if (mode == 0)
	//{
	//    for (double b : beta) std::cout << b << std::endl;
	//}
	//else
	//{
	//    for (double x : betaPJ) std::cout << x << std::endl;
	//}

    return 0;
}