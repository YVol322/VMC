#pragma once
#include <memory>
#include <cmath>
#include <vector>

#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

class Sampler {
public:
    Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength,
        unsigned int numberOfMetropolisSteps);


    void sample(bool acceptedStep, class System* system);
    void printOutputToTerminal(class System& system);
    void computeAverages();
    void setEnergy(double en);
    void setTime(double t);
    double getEnergy() { return m_energy; }
    std::vector<double> getEnergyrij() { return m_Elrij; }
    std::vector<double> getrij() { return m_rij; }

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;
    double m_energy = 0;
    double m_cumulativeEnergy = 0;
    double m_stepLength = 0;
    double m_time = 0;
    std::vector<double> m_rij;
    std::vector<double> m_Elrij;
    std::vector<double> m_cumulativeEnergyrij;
    std::vector<double> m_cumulativerij;
};