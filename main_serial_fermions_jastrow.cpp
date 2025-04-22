#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
#include "WaveFunctions/fermionsjastrow2autodiff.h"
#include "WaveFunctions/fermionsjastrow2.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "Solvers/metropolis.h"
#include "Math/random.h"
#include "particle.h"
#include "sampler.h"
#include <autodiff/forward/dual.hpp>


using namespace std;
using namespace std::chrono;


int main() {
    
    int seed = 2025;

    unsigned int numberOfDimensions = 2;
    unsigned int numberOfParticles = 2;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e4;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e3;

    double omega = 1.0;
    double alpha = 0.5;
    double beta = 0;

    double stepLength = 1;

	double learning_rate = 1e-1;
	double stop_at = 1e-2;
	double max_iters = 1000;
	double iter = 0;
	double grad_beta = 1;
	
	while(iter < max_iters && abs(grad_beta) > stop_at)
	{
		auto rng = std::make_unique<Random>(seed);
		auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);
		auto system = std::make_unique<System>(
        std::make_unique<HarmonicOscillator>(omega),
        std::make_unique<FermionsJastrow2>(alpha, beta),
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

    	double duration = duration_cast<seconds>(stop - start).count();

		sampler -> setTime(duration);
    	//sampler -> printOutputToTerminal(*system);

		double mean_rij = (sampler -> getrij()).at(0);
		double mean_El = sampler -> getEnergy();
		double mean_El_times_rij = (sampler -> getEnergyrij()).at(0);

		grad_beta = 2 * (mean_El_times_rij - mean_El * mean_rij);

		beta -= learning_rate * grad_beta;

		std::cout << "Iteration " << iter
          << ", beta = " << beta
          << ", energy = " << mean_El
          << ", grad = " << grad_beta << std::endl;


		iter++;
	}

    return 0;
}