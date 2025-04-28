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
    m_energy = 0;
    m_cumulativeEnergy = 0;
    m_stepLength = stepLength;
    m_numberOfAcceptedSteps = 0;
    
    m_nPairs = m_numberOfParticles * (m_numberOfParticles - 1) / 2;

    m_rij = std::vector<double>(m_nPairs, 0.0);
    m_Elrij = std::vector<double>(m_nPairs, 0.0);
    m_cumulativeEnergyrij = std::vector<double>(m_nPairs, 0.0);
    m_cumulativerij = std::vector<double>(m_nPairs, 0.0);

    m_O = std::vector<double>(1, 0.0);
    m_ElO = std::vector<double>(1, 0.0);
    m_cumulativeElO = std::vector<double>(1, 0.0);
    m_cumulativeO = std::vector<double>(1, 0.0);
}


void Sampler::sample(bool acceptedStep, System* system) {
    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;
    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;

    int N = m_numberOfParticles;
    int pairIndex = 0;

    double O = 0.0;
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            double rij = system->computerij(i, j);
            m_cumulativerij[pairIndex] += rij;
            m_cumulativeEnergyrij[pairIndex] += rij * localEnergy;

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
    m_cumulativeO[0] += O;
    m_cumulativeElO[0] += O * localEnergy;
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
    cout << " Numberical Energy : " << m_energy << endl;
    cout << endl;

    cout << " Analytical Energy for bosons : " << 0.5 * m_numberOfDimensions * m_numberOfParticles * (pa.back() + 1/(4* pa.back())) << endl;
    cout << " Analytical Energy for fermions : " << 10 << endl;
    cout << endl;

    cout << " Algo runtime: " << m_time << " seconds" << endl;
}

void Sampler::computeAverages() {
    m_energy = m_cumulativeEnergy / m_numberOfMetropolisSteps;

    for (int i = 0; i < m_nPairs; i++)
    {
        m_rij[i] = m_cumulativerij[i] / m_numberOfMetropolisSteps;
        m_Elrij[i] = m_cumulativeEnergyrij[i] / m_numberOfMetropolisSteps;
    }
    m_ElO[0] = m_cumulativeElO[0] / m_numberOfMetropolisSteps;
    m_O[0] = m_cumulativeO[0] / m_numberOfMetropolisSteps;
}

void Sampler::setEnergy(double en) {
    m_energy = en;
}

void Sampler::setTime(double t) {
    m_time = t;
}