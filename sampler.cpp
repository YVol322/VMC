#include <iostream>         // Include the C++ input-output stream library.

#include "sampler.h"        // Include "sampler" header file with declarations.

using std::cout;
using std::endl;


// Constructor for the Sampler class. Initializes all necessary parameters for the sampling process 
// in a Variational Monte Carlo (VMC) simulation. This constructor sets up the number of particles, 
// dimensions, step length, and the number of Metropolis steps, along with initializing variables 
// to store the cumulative values of various observables for optimization purposes.
//
// Input:   unsigned int numberOfParticles - the number of particles in the system;
//          unsigned int numberOfDimensions - the number of dimensions in the system;
//          double stepLength - the step length used in the Metropolis algorithm;
//          unsigned int numberOfMetropolisSteps - the total number of Metropolis steps in the simulation.
Sampler::Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength,
        unsigned int numberOfMetropolisSteps)
{
    m_stepNumber = 0;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    m_stepLength = stepLength;
    m_numberOfAcceptedSteps = 0;
    
    m_nPairs = m_numberOfParticles * (m_numberOfParticles - 1) / 2;

    m_energy = 0;
    m_cumulativeEnergy = 0;

    m_O1alpha = 0;
    m_cumulativeO1alpha = 0;

    m_O2alpha = 0;
    m_cumulativeO2alpha = 0;

    m_O1Jastrow = std::vector<double>(m_nPairs, 0.0);
    m_cumulativeO1Jastrow = std::vector<double>(m_nPairs, 0.0);
    
    m_O2Jastrow = std::vector<double>(m_nPairs, 0.0);
    m_cumulativeO2Jastrow = std::vector<double>(m_nPairs, 0.0);

    m_O1Pade = 0;
    m_cumulativeO1Pade = 0;

    m_O2Pade = 0;
    m_cumulativeO2Pade = 0;
}

// Samples all observables, including the local energy, O1, and O2 for alpha, beta, and beta_ij optimization.
//
// Input:   bool acceptedStep - indicates if the step was accepted in the Metropolis algorithm;
//          class System* system - pointer to the system object that contains the simulation state.
void Sampler::sample(bool acceptedStep, System* system)
{
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;     // Sample local energy.

    double r2 = system -> computer2();
    m_cumulativeO1alpha -= r2;                  // Sample O1alpha for alpha optimization.
    m_cumulativeO2alpha -= r2 * localEnergy;    // Sample O2alpha for alpha optimization.

    int N = m_numberOfParticles;
    int pairIndex = 0;

    double O = 0.0;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            double rij = system -> computerij(i, j);
            m_cumulativeO1Jastrow[pairIndex] += rij;                // Sample O1Jastrow for beta_ij optimization.
            m_cumulativeO2Jastrow[pairIndex] += rij * localEnergy;  // Sample O2Jastrow for beta_ij optimization.

            double aij;
            if ((i < N/2 && j >= N/2) || (i >= N/2 && j < N/2))
            {
                aij = 1.0;
            }
            else
            {
                aij = 1.0/3.0;
            }

            double beta = system -> getWaveFunctionParameters().at(0);
            double t = 1.0 + beta * rij;

            O += -aij * rij * rij / (t * t);

            pairIndex++;
        }
    }
    m_cumulativeO1Pade += O;                    // Sample O1Pade for beta optimization.
    m_cumulativeO2Pade += O * localEnergy;      // sample O2Pade for beta optimization.

    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;
}


// Prints the output to the terminal, typically the results of the VMC simulation.
//
// Input:   class System& system - reference to the system object to extract the necessary information.
void Sampler::printOutputToTerminal(System& system) {
    auto pa = system.getWaveFunctionParameters();
    auto p = pa.size();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << m_numberOfParticles << endl;
    cout << " Number of dimensions : " << m_numberOfDimensions << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(m_numberOfMetropolisSteps) << endl;
    cout << " Step length used : " << m_stepLength << endl;
    cout << " Ratio of accepted steps: " << ((double) m_numberOfAcceptedSteps) / ((double) m_numberOfMetropolisSteps) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (unsigned int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Numerical Energy : " << m_energy << endl;
    cout << endl;

    cout << " Analytical Energy for bosons : " << 0.5 * m_numberOfDimensions * m_numberOfParticles * (pa.back() + 1/(4* pa.back())) << endl;
    cout << endl;

    cout << " Algo runtime: " << m_time << " seconds" << endl;
}


// Computes the averages of all sampled variables after the sampling is complete.
void Sampler::computeAverages()
{
    m_energy = m_cumulativeEnergy / m_numberOfMetropolisSteps;

    m_O1alpha = m_cumulativeO1alpha / m_numberOfMetropolisSteps;
    m_O2alpha = m_cumulativeO2alpha / m_numberOfMetropolisSteps;

    for (int i = 0; i < m_nPairs; i++)
    {
        m_O1Jastrow[i] = m_cumulativeO1Jastrow[i] / m_numberOfMetropolisSteps;
        m_O2Jastrow[i] = m_cumulativeO2Jastrow[i] / m_numberOfMetropolisSteps;
    }
    m_O1Pade = m_cumulativeO1Pade / m_numberOfMetropolisSteps;
    m_O2Pade = m_cumulativeO2Pade / m_numberOfMetropolisSteps;
}


// Sets the VMC energy to be used in output (necessary for parallel algorithms).
//
// Input:   double en - the value of the energy to set.
void Sampler::setEnergy(double en)
{
    m_energy = en;
}


// Sets the VMC time to be used in output (necessary for parallel algorithms).
//
// Input:   double t - the runtime of the VMC simulation.
void Sampler::setTime(double t)
{
    m_time = t;
}