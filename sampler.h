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
    double getO1alpha() { return m_O1alpha; }
    double getO2alpha() { return m_O2alpha; }
    std::vector<double> getO1Jastow() { return m_O1Jastrow; }
    std::vector<double> getO2Jastow() { return m_O2Jastrow; }
    double getO1Pade() { return m_O1Pade; }
    double getO2Pade() { return m_O2Pade; }

private:
    unsigned int m_stepNumber = 0;
    unsigned int m_numberOfMetropolisSteps = 0;
    unsigned int m_numberOfParticles = 0;
    unsigned int m_numberOfDimensions = 0;
    unsigned int m_numberOfAcceptedSteps = 0;

    double m_stepLength = 0;
    double m_time = 0;
    int m_nPairs = 0;

    double m_energy = 0;
    double m_cumulativeEnergy = 0;

    double m_O1alpha = 0;
    double m_cumulativeO1alpha = 0;

    double m_O2alpha = 0;
    double m_cumulativeO2alpha = 0;


    std::vector<double> m_O1Jastrow;
    std::vector<double> m_cumulativeO1Jastrow;

    std::vector<double> m_O2Jastrow;
    std::vector<double> m_cumulativeO2Jastrow;

    double m_O1Pade;
    double m_cumulativeO1Pade = 0;

    double m_O2Pade;
    double m_cumulativeO2Pade = 0;
};