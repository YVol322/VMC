#pragma once

#include <cassert>              // Include the C++ assert library for runtime assertions
#include <memory>               // Include the C++ memory library for using std::unique_ptr
#include <vector>               // Include the C++ vector library for using std::vector<>

#include "../system.h"          // Include "system" header file with declarations.
#include "../particle.h"        // Include "particle" header file with declarations.


// Declaration of the WaveFunction class. This class is designed to compute a trial wave function,
// the sum of second derivatives, and the quantum force.
// Each subclass of this class implements its own ansatz. 
// This class defines the functions that can be used by its subclasses.
class WaveFunction
{
    public:
        // Destructor of the WaveFunction class. Ensures proper cleanup of derived class objects.
        virtual ~WaveFunction() = default;
    

        // Helper function to get the number of parameters for the wave function.
        //
        // Output: int - the number of parameters.
        int getNumberOfParameters() { return m_numberOfParameters; }


        // Helper function to get the vector of parameters for the wave function.
        //
        // Output: const std::vector<double>& - reference to the vector of parameters.
        const std::vector<double>& getParameters() { return m_parameters; }


        // Computes the sum of squared coordinates for the particle specified by part_index.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int part_index - index of the particle for which to compute r^2.
        //
        // Output:  double - the sum of squared coordinates for the specified particle.
        double r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index);


        // Computes the relative distance between particles i and j.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int i - index of the first particle;
        //          int j - index of the second particle.
        //
        // Output:  double - the relative distance between particles i and j.
        double r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j);


        // Transforms a 2D index (i, j) into a 1D index k.
        //
        // Input:   int i - index of the first particle;
        //          int j - index of the second particle.
        //
        // Output:  int - the corresponding 1D index.
        int BetaIndex(int i, int j);


        // Computes the constant a (for Pade-Jastrow ansatz) for particles i and j.
        //
        // Input:   int i - index of the first particle;
        //          int j - index of the second particle.
        //
        // Output:  double - the constant a for particles i and j.
        double a_ij(int i, int j);


        // Computes the Green's function ratio between two particle configurations.
        //
        // Input:   std::vector<double>& Rnew - new configuration of particle positions;
        //          std::vector<double>& Rold - old configuration of particle positions;
        //          double dt - time step used for the update;
        //          std::vector<double>& Fold - forces for the old configuration;
        //          std::vector<double>& Fnew - forces for the new configuration.
        //
        // Output:  double - the computed Green's function ratio.
        double GreensFunctionRatio(std::vector<double>& Rnew, std::vector<double>& Rold, 
                                   double dt, std::vector<double>& Fold, std::vector<double>& Fnew);


        // Computes the ansatz for the wave function. This function is virtual and must be implemented by subclasses.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed wave function value.
        virtual double evaluate(std::vector<std::unique_ptr<class Particle>>& particles) = 0;


        // Computes the sum of second derivatives of the wave function. This function is virtual and must be implemented by subclasses.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
        //
        // Output:  double - the computed sum of second derivatives.
        virtual double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles) = 0;


        // Computes the quantum force for the specified particle. This function is virtual and must be implemented by subclasses.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int i - index of the particle for which to compute the quantum force.
        //
        // Output:  std::vector<double> - the computed quantum force for particle i.
        virtual std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int i) = 0;
    
    protected:
        int m_numberOfParameters = 0;         // Number of parameters for the wave function.
        std::vector<double> m_parameters = std::vector<double>(); // Vector of parameters for the wave function.
        int m_particles = 0;                  // Number of particles in the system.
        int m_mode = 0;                       // Mode for the ansatz: 0 - Jastrow ansatz, 1 - Pade-Jastrow.
        double m_omega = 0;                   // Angular frequency for the harmonic oscillator (if applicable).
        double m_sqrt_om = 0;                 // Square root of the angular frequency.
};