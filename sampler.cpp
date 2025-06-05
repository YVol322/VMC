#include <iostream>

#include "sampler.h"

using std::cout;
using std::endl;


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


void Sampler::sample(bool acceptedStep, System* system)
{
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;

    double r2 = system -> computer2();
    m_cumulativeO1alpha -= r2;
    m_cumulativeO2alpha -= r2 * localEnergy;

    int N = m_numberOfParticles;
    int pairIndex = 0;

    double O = 0.0;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            double rij = system -> computerij(i, j);
            m_cumulativeO1Jastrow[pairIndex] += rij;
            m_cumulativeO2Jastrow[pairIndex] += rij * localEnergy;

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
    m_cumulativeO1Pade += O;
    m_cumulativeO2Pade += O * localEnergy;

    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;
}


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
    cout << " Analytical Energy for fermions : " << 10 << endl;
    cout << endl;

    cout << " Algo runtime: " << m_time << " seconds" << endl;
}


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

void Sampler::setEnergy(double en)
{
    m_energy = en;
}

void Sampler::setTime(double t)
{
    m_time = t;
}