#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(double alpha);
    double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);
    double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
};
