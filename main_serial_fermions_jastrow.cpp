#include <iostream>
#include <vector>
#include <memory>
#include <chrono>

#include "system.h"
#include "WaveFunctions/fermionsjastrownumerical.h"
#include "WaveFunctions/fermionsjastrow.h"
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
    unsigned int numberOfParticles = 12;
    unsigned int numberOfMetropolisSteps = (unsigned int) 1e3;
    unsigned int numberOfEquilibrationSteps = (unsigned int) 1e2;

    double omega = 1.0;
    double alpha = 0.5;

	int nPairs = numberOfParticles * (numberOfParticles - 1) / 2;
    std::vector<double> beta(nPairs, 0);
    //std::vector<double> beta(nPairs, 0.145819);

	//beta[0] = 0.142175;
	//beta[1] = 0.145953;
	//beta[2] = 0.19525;
	//beta[3] = 0.20664;
	//beta[4] = 0.211661;
	//beta[5] = 0.127557;
	//beta[6] = 0.224082;
	//beta[7] = 0.162106;
	//beta[8] = 0.223895;
	//beta[9] = 0.157206;
	//beta[10] = 0.205652;
	//beta[11] = 0.209072;
	//beta[12] = 0.133737;
	//beta[13] = 0.151694;
	//beta[14] = 0.132747;


    double stepLength = 1;

	double learning_rate = 1e-2;
	double stop_at = 5e-2;
	double max_iters = 1000;
	double iter = 0;
    std::vector<double> grad_beta(nPairs, 1);
	double l2_norm = 1.0;
	
	while(iter < max_iters && l2_norm > stop_at)
	{
		auto rng = std::make_unique<Random>(seed);
		auto particles = setupRandomUniformInitialState(numberOfDimensions, numberOfParticles, *rng);
		auto system = std::make_unique<System>(
        std::make_unique<HarmonicOscillator>(omega),
        std::make_unique<FermionsJastrow>(alpha, beta),
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
	//for(int i = 0; i < nPairs; i++) std::cout << beta[i] << std::endl;

    return 0;
}