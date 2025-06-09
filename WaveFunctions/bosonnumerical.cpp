#include "bosonnumerical.h"     // Include "bosonnumerical" header file with declarations.


// Constructor of the BosonNumerical class. Initializes the parameters for the boson system.
// It asserts that the variational parameter alpha is non-negative and initializes the necessary parameters
// including the number of particles, omega (angular frequency), and its square root.
// 
// Input:   double alpha - variational parameter alpha;
//          int n_particles - number of particles in the system;
//          double omega - angular frequency of the system.
BosonNumerical::BosonNumerical(double alpha, int n_particles, double omega)
{
    assert(alpha >= 0);  // Ensure the variational parameter alpha is non-negative.

    m_numberOfParameters = 1;  // Only one parameter, alpha.

    m_parameters.reserve(m_numberOfParameters);  // Reserve space for the parameter.
    m_parameters.push_back(alpha);  // Store alpha as the variational parameter.

    m_particles = n_particles;  // Set the number of particles.

    m_omega = omega;  // Set the angular frequency.
    m_sqrt_om = sqrt(omega);  // Calculate the square root of omega.
}


// BosonNumerical subclass specific function that fills the state vector x with the particle positions.
// This function is used to generate the vector of particle positions required for further calculations.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
// 
// Output:  VectorXvar - the state vector x, containing the positions of all particles.
VectorXvar BosonNumerical::fill_x(std::vector<std::unique_ptr<class Particle>>& particles)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    VectorXvar x(n_particles * n_dimensions);  // Initialize the state vector with size equal to the total number of coordinates.
    
    // Fill the state vector x with the positions of the particles.
    for(int i = 0; i < n_particles; i++)
    {
        std::vector<double> ri = particles[i] -> getPosition();
        for(int j = 0; j < n_dimensions; j++)
        {
            x(i * n_dimensions + j) = ri[j];
        }
    }
    return x;
}


// BosonNumerical subclass specific function that computes the trial wave function (Psi_T) 
// using the state vector x (which contains the particle positions).
//
// Input:   VectorXvar& x - the state vector containing the particle positions.
// 
// Output:  var - the value of the trial wave function (Psi_T) for the given particle configuration.
var BosonNumerical::PsiT(VectorXvar& x)
{
    int n_dimensions = x.size() / m_particles;

    var sum = 0;
    
    // Sum the squared coordinates of each particle to calculate r^2.
    for(int i = 0; i < m_particles; i++)
    {
        for(int j = 0; j < n_dimensions; j++)
        {
            var pos = x(i * n_dimensions + j);
            sum += pos * pos;
        }
    }
    
    var alpha = m_parameters.back();

    // Return the trial wave function: Psi_T = exp(-alpha * omega * r^2). The output type is 'var',
    // allowing the autodiff library to compute derivatives with respect to the entries of the vector x.
    return exp(-alpha * m_omega * sum);
}


// BosonNumerical subclass specific function that computes the Laplacian of Psi_T over Psi_T 
// using the state vector x and automatic differentiation.
//
// Input:   VectorXvar& x - the state vector containing the particle positions.
// 
// Output:  double - the computed Laplacian of Psi_T divided by Psi_T for all particles in the system.
double BosonNumerical::LaplPsiTOverPsiT(VectorXvar& x)
{
    var psiT = PsiT(x);  // Compute Psi_T using the state vector x.

    Eigen::VectorXd g;  // Gradient vector (not used here but required for Hessian calculation).
    Eigen::MatrixXd H = hessian(psiT, x, g);  // Compute the Hessian matrix of Psi_T with respect to the state vector x.
    
    Eigen::VectorXd LaplVector = H.diagonal();  // Extract the diagonal of the Hessian, which corresponds to the Laplacian.

    double sum = val(LaplVector.sum());  // Sum the diagonal elements and extract the value.

    return sum / val(psiT);  // Return the Laplacian of Psi_T divided by Psi_T.
}


// WaveFunction class function that evaluates the wave function for the given particle configuration 
// using the bosonic ansatz and automatic differentiation. It calls the PsiT function for evaluation.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
// 
// Output:  double - the value of the wave function for the given particle configuration using the bosonic ansatz.
double BosonNumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    VectorXvar x = fill_x(particles);  // Fill the state vector x with particle positions.

    return val(PsiT(x));
}


// WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
// for the bosonic ansatz, which is required for calculating the local energy using automatic differentiation.
// It calls the LaplPsiTOverPsiT function for the Laplacian.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
// 
// Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
double BosonNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);  // Fill the state vector x with particle positions.

    return LaplPsiTOverPsiT(x);
}


// WaveFunction class function that computes the quantum force for a specific particle in the system
// using the bosonic ansatz and automatic differentiation. The quantum force is used in the Metropolis-Hastings sampling algorithm.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int n - index of the particle for which to compute the quantum force.
// 
// Output:  std::vector<double> - the quantum force for the specified particle.
std::vector<double> BosonNumerical::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    VectorXvar x = fill_x(particles);  // Fill the state vector x with particle positions.
    var psiT = PsiT(x);

    int d = x.size() / m_particles;
    std::vector<double> F(d);
    
    // Compute the quantum force for the specified particle in all dimensions.
    for(int i = 0; i < d; i++)
    {
        auto [gradx] = derivatives(psiT, wrt(x(d * n + i)));  // Compute the derivative of Psi_T with respect to the position.
        F[i] = 2 * val(gradx) / (val(psiT));
    }
    
    return F;
}
