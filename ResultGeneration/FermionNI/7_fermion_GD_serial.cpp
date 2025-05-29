#include <iostream>
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

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 2;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e3;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e2;

    double omega = 1.0;
    double alpha = 0.6;

	double stepLength = 1;
	double learning_rate = 1e-1;
	double stop_at = 1e-3;
	double max_iters = 1000;
	double iter = 0;
	double l2_norm = 4.1;
	
	while(iter < max_iters && l2_norm > stop_at)
	{
		auto rng = std::make_unique<Random>(seed);
		auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);
		auto system = std::make_unique<System>(
        std::make_unique<HarmonicOscillator>(omega, 0),
		std::make_unique<FermionNI>(alpha, numberOfParticles),
        std::make_unique<Metropolis>(std::move(rng)),
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

		double El = sampler -> getEnergy();
		double O1 = sampler -> getO1alpha();
		double O2 = sampler -> getO2alpha();
	
		double grad = 2 * (O2 - El * O1);
	
		alpha += learning_rate * grad;
		l2_norm = abs(grad);
		std::cout << "iter = " << iter << std::endl;
        std::cout << "local en = "<< El << std::endl;
        std::cout << "grad = "<< grad << std::endl;
        std::cout << "alpha = "<< alpha << std::endl;
        std::cout << "time = "<< duration.count() << std::endl;

		iter++;
	}

    return 0;
}