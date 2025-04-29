#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
#include "WaveFunctions/fermionsjastrownumerical.h"
#include "WaveFunctions/fermionsjastrow.h"
#include "WaveFunctions/boson.h"
#include "WaveFunctions/bosonnumerical.h"
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
    unsigned int numberOfParticles = 6;
	int mode = 1;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e3;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e2;

    double omega = 1.0;
    double alpha = 0.5;

	int nPairs = numberOfParticles * (numberOfParticles - 1) / 2;
    std::vector<double> beta(nPairs, 0.3);
	std::vector<double> betaPJ(1, 0.446563);


    double stepLength = 1;

	double learning_rate = 1e-2;
	double stop_at = 1e-2;
	double max_iters = 1000;
	double iter = 0;
    std::vector<double> grad_beta(nPairs, 1);
	double l2_norm = 4.1;

    std::vector<double> grad_betaPJ(1, 1);
	
	while(iter < max_iters && l2_norm > stop_at)
	{
		auto rng = std::make_unique<Random>(seed);
		auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);
		auto system = std::make_unique<System>(
        std::make_unique<HarmonicOscillator>(omega),
        //std::make_unique<FermionsJastrow>(alpha, (mode == 0 ? beta : betaPJ), mode, numberOfParticles),
		std::make_unique<BosonNumerical>(alpha, betaPJ, mode, numberOfParticles),
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
		if(mode == 0)
		{
			double mean_rij;
			double mean_El;
			double mean_El_times_rij;
			for (int i = 0; i < nPairs; i++)
			{
				mean_rij = (sampler -> getrij()).at(i);
				mean_El = sampler -> getEnergy();
				mean_El_times_rij = (sampler -> getEnergyrij()).at(i);

				grad_beta.at(i) = 2 * (mean_El_times_rij - mean_El * mean_rij);

				beta.at(i) -= learning_rate * grad_beta.at(i);
			}

			l2_norm = 0.0;
			for (double g : grad_beta)
			{
			    l2_norm += g * g;
			}
			l2_norm = sqrt(l2_norm);

			std::cout << "Iteration " << iter
        	  << ", beta = " << beta[0]
        	  << ", energy = " << mean_El
        	  << ", grad = " << l2_norm << std::endl;

			iter++;
		}
		else
		{
			double mean_O;
			double mean_El_times_O;
			double mean_El;
			
			mean_O = (sampler -> getO()).at(0);
			mean_El = sampler -> getEnergy();
			mean_El_times_O = (sampler -> getEnergyO()).at(0);
	
			grad_betaPJ.at(0) = 2 * (mean_El_times_O - mean_El * mean_O);
	
			betaPJ.at(0) -= learning_rate * grad_betaPJ.at(0);
	
			l2_norm = 0.0;
			for (double g : grad_betaPJ)
			{
			    l2_norm += g * g;
			}
			l2_norm = sqrt(l2_norm);

			std::cout << "Iteration " << iter
        	  << ", beta = " << betaPJ[0]
        	  << ", energy = " << mean_El
        	  << ", grad = " << l2_norm << std::endl;

			iter++;
		}
	}

	if (mode == 0)
	{
	    for (double b : beta) std::cout << b << std::endl;
	}
	else
	{
	    for (double x : betaPJ) std::cout << x << std::endl;
	}

    return 0;
}