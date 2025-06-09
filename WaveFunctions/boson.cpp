#include "boson.h"      // Include "boson" header file with declarations.


// Constructor of the Boson class. Initializes the parameters for the boson system.
// It asserts that the variational parameter alpha is non-negative, and initializes the necessary parameters.
// 
// Input:   double alpha - variational parameter alpha;
//          int n_particles - number of particles in the system;
//          double omega - angular frequency of the system.
Boson::Boson(double alpha, int n_particles, double omega)
{
    assert(alpha >= 0);  // Ensure the variational parameter alpha is non-negative.

    m_numberOfParameters = 1;  // Only one parameter, alpha.

    m_parameters.reserve(m_numberOfParameters);  // Reserve space for the parameter.
    m_parameters.push_back(alpha);  // Store alpha as the variational parameter.

    m_particles = n_particles;  // Set the number of particles.

    m_omega = omega;  // Set the angular frequency.
    m_sqrt_om = sqrt(omega);  // Calculate the square root of omega.
}


// Boson subclass specific function that computes the trial wave function (Psi_T).
// The wave function is given by Psi_T = exp(-alpha * omega * r^2), where r^2 is the sum of squared coordinates of all particles.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the value of the trial wave function (Psi_T) for the given particle configuration.
double Boson::PsiT(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();  // Get the variational parameter alpha.
    double argument = 0;

    // Sum the squared coordinates (r^2) of all particles.
    for(int i = 0; i < m_particles; i++)
    {
        argument += r_squared(particles, i);
    }

    return exp(-alpha * m_omega * argument);
}


// Boson subclass specific function that computes the sum of Laplacians of Psi_T divided by Psi_T for all particles.
// This is used in the local energy computation.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
double Boson::LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    double argument = 0;

    Particle& particle_i = *(particles[0]);
    int n_dimensions = particle_i.getNumberOfDimensions();

    // Sum the squared coordinates (r^2) of all particles.
    for(int i = 0; i < m_particles; i++)
    {
        argument += r_squared(particles, i);
    }

    return (-2 * n_dimensions * m_particles * alpha * m_omega + 4 * alpha * alpha * m_omega * m_omega * argument);
}


// WaveFunction class function that evaluates the wave function for the given particle configuration
// using the bosonic ansatz. It calls PsiT for evaluation.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the value of the wave function for the given particle configuration using the bosonic ansatz.
double Boson::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    return PsiT(particles);
}


// WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
// for the bosonic ansatz, which is required for calculating the local energy. It calls LaplPsiTOverPsiT for the Laplacian.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
double Boson::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    return LaplPsiTOverPsiT(particles);  // Return the computed sum of second derivatives (Laplacian).
}


// WaveFunction class function that computes the quantum force for a specific particle in the system
// using the bosonic ansatz. The quantum force is used in the Metropolis-Hastings sampling algorithm.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int n - index of the particle for which to compute the quantum force.
//
// Output:  std::vector<double> - the quantum force for the specified particle.
std::vector<double> Boson::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    double alpha = m_parameters.back();
    double d = particles[n]->getNumberOfDimensions();
    std::vector<double> r = particles[n]->getPosition();
    std::vector<double> F(d);  // Initialize the force vector.

    // Compute the quantum force for the specified particle in all dimensions.
    for(int i = 0; i < d; i++)
    {
        F[i] = -4 * alpha * m_omega * r[i];
    }

    return F;
}