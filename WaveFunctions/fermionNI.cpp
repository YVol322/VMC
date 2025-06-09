#include "fermionNI.h"      // Include "bosonnumerical" header file with declarations.


// Constructor of the FermionNI class. Initializes the parameters for the fermion system.
// It sets the number of parameters, stores the variational parameter alpha, 
// and initializes the number of particles and angular frequency (omega).
//
// Input:   double alpha - variational parameter alpha;
//          int n_particles - number of particles in the system;
//          double omega - angular frequency of the system.
FermionNI::FermionNI(double alpha, int n_particles, double omega)
{
    assert(alpha >= 0);  // Ensure the variational parameter alpha is non-negative.

    m_numberOfParameters = 1;  // Only one parameter, alpha.

    m_parameters.reserve(m_numberOfParameters);  // Reserve space for the parameter.
    m_parameters.push_back(alpha);  // Store alpha as the variational parameter.

    m_particles = n_particles;  // Set the number of particles.

    m_omega = omega;  // Set the angular frequency.
    m_sqrt_om = sqrt(omega);  // Calculate the square root of omega.
}


// Computes the value of Psi1 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi1.
//
// Output:  double - the value of Psi1 for the given particle.
double FermionNI::Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    double r2 = r_squared(particles, part_idx);

    return exp(-alpha * m_omega * r2);
}


// Computes the value of Psi2 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi2.
//
// Output:  double - the value of Psi2 for the given particle.
double FermionNI::Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double x = particles[part_idx] -> getPosition()[0];

    return x * Psi1(particles, part_idx);
}


// Computes the value of Psi3 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi3.
//
// Output:  double - the value of Psi3 for the given particle.
double FermionNI::Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double y = particles[part_idx] -> getPosition()[1];

    return y * Psi1(particles, part_idx);
}


// Computes the value of Psi4 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi4.
//
// Output:  double - the value of Psi4 for the given particle.
double FermionNI::Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return x * y * Psi1(particles, part_idx);
}


// Computes the value of Psi5 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi5.
//
// Output:  double - the value of Psi5 for the given particle.
double FermionNI::Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double x = particles[part_idx] -> getPosition()[0];

    return (x * x - 1) * Psi1(particles, part_idx);
}


// Computes the value of Psi6 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi6.
//
// Output:  double - the value of Psi6 for the given particle.
double FermionNI::Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double y = particles[part_idx] -> getPosition()[1];

    return (y * y - 1) * Psi1(particles, part_idx);
}


// Computes the gradient of Psi1 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi1.
//
// Output:  std::vector<double> - the gradient of Psi1 for the given particle.
std::vector<double> FermionNI::GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi1 = {-2 * alpha * m_omega * x * psi1, -2 * alpha * m_omega * y * psi1};

    return grad_psi1;
}


// Computes the gradient of Psi2 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi2.
//
// Output:  std::vector<double> - the gradient of Psi2 for the given particle.
std::vector<double> FermionNI::GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi2 = {(1 - 2 * alpha * m_omega * x * x) * psi1, -2 * alpha * m_omega * y * x * psi1};

    return grad_psi2;
}


// Computes the gradient of Psi3 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi3.
//
// Output:  std::vector<double> - the gradient of Psi3 for the given particle.
std::vector<double> FermionNI::GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi3 = {-2 * alpha * m_omega * x * y * psi1, (1 - 2 * alpha * m_omega * y * y) * psi1};

    return grad_psi3;
}


// Computes the gradient of Psi4 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi4.
//
// Output:  std::vector<double> - the gradient of Psi4 for the given particle.
std::vector<double> FermionNI::GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi4 = {y * psi1 * (1 - 2 * alpha * m_omega * x * x), x * psi1 * (1 - 2 * alpha * m_omega * y * y)};

    return grad_psi4;
}


// Computes the gradient of Psi5 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi5.
//
// Output:  std::vector<double> - the gradient of Psi5 for the given particle.
std::vector<double> FermionNI::GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi5 = {2 * x * psi1 * (1 - alpha * m_omega * x * x + alpha * m_omega), -2 * alpha * m_omega * y * (x * x - 1) * psi1};

    return grad_psi5;
}


// Computes the gradient of Psi6 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi6.
//
// Output:  std::vector<double> - the gradient of Psi6 for the given particle.
std::vector<double> FermionNI::GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi6 = {-2 * alpha * m_omega * x * (y * y - 1) * psi1, 2 * y * psi1 * (1 - alpha * m_omega * y * y + alpha * m_omega)};

    return grad_psi6;
}


// Computes the Laplacian of Psi1 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi1.
//
// Output:  double - the Laplacian of Psi1 for the given particle.
double FermionNI::LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-4 * alpha * m_omega+ 4 * alpha * alpha * m_omega * m_omega * (x * x + y * y) ) * Psi1(particles, part_idx);
}


// Computes the Laplacian of Psi2 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi2.
//
// Output:  double - the Laplacian of Psi2 for the given particle.
double FermionNI::LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-8 * alpha * m_omega * x + 4 * alpha * alpha * m_omega * m_omega * x * (x * x + y * y) ) * Psi1(particles, part_idx);
}


// Computes the Laplacian of Psi3 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi3.
//
// Output:  double - the Laplacian of Psi3 for the given particle.
double FermionNI::LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-8 * alpha * m_omega * y + 4 * alpha * alpha * m_omega * m_omega * y * (x * x + y * y) ) * Psi1(particles, part_idx);
}


// Computes the Laplacian of Psi4 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi4.
//
// Output:  double - the Laplacian of Psi4 for the given particle.
double FermionNI::LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return x * y * Psi1(particles, part_idx) * (4 * alpha * alpha * m_omega * m_omega *(x * x + y * y)
                                                - 12 * alpha * m_omega);
}


// Computes the Laplacian of Psi5 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi5.
//
// Output:  double - the Laplacian of Psi5 for the given particle.
double FermionNI::LapliPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    double x2 = x * x;
    double y2 = y * y;
    double common = x2 - 1;

    double laplacian_psi5 = 2 * psi1 *
    (-4 * alpha * m_omega * x2 + alpha * m_omega * common * (2 * alpha * m_omega * x2 - 1)
    + alpha * m_omega * common * (2 * alpha * m_omega * y2 - 1) + 1);

    return laplacian_psi5;
}


// Computes the Laplacian of Psi6 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi6.
//
// Output:  double - the Laplacian of Psi6 for the given particle.
double FermionNI::LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    double x2 = x * x;
    double y2 = y * y;
    double common = y2 - 1;

    double laplacian_psi6 = 2 * psi1 *
    (-4 * alpha * m_omega * y2 + alpha * m_omega * common * (2 * alpha * m_omega * x2 - 1)
    + alpha * m_omega * common * (2 * alpha * m_omega * y2 - 1) + 1);

    return laplacian_psi6;
}


// Computes the spin-up Slater determinant, its gradient component, or Laplacian depending on the arguments.
        //
        // Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
        //          int particles_set - number of sets of particles: 0 - spin-up, 1 - spin-down;
        //          int row_changed - row index of the particle that was changed;
        //          int der_order - order of the derivatives to compute;
        //          int grad_comp - the gradient component to compute (if applicable).
        //
        // Output:  double - the computed value based on the requested derivative order (Slater determinant, gradient, or Laplacian).
double FermionNI::SD
(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp)
{
    int size = m_particles / 2;
    double det = 0;

    if(m_particles == 6)
    {
        Eigen::MatrixXd A(size, size);
        for(int i = 0; i < size; i++)
        {
            if (i == row_changed)
            {
                if (der_order == 0)
                {
                    A(i, 0) = Psi1(particles, i + size * particles_set);
                    A(i, 1) = Psi2(particles, i + size * particles_set);
                    A(i, 2) = Psi3(particles, i + size * particles_set);
                }
                else if (der_order == 1)
                {
                    A(i, 0) = GradiPsi1(particles, i + size * particles_set).at(grad_comp);
                    A(i, 1) = GradiPsi2(particles, i + size * particles_set).at(grad_comp);
                    A(i, 2) = GradiPsi3(particles, i + size * particles_set).at(grad_comp);
                }
                else if (der_order == 2)
                {
                    A(i, 0) = LapliPsi1(particles, i + size * particles_set);
                    A(i, 1) = LapliPsi2(particles, i + size * particles_set);
                    A(i, 2) = LapliPsi3(particles, i + size * particles_set);
                }
            }
            else
            {
                A(i, 0) = Psi1(particles, i + size * particles_set);
                A(i, 1) = Psi2(particles, i + size * particles_set);
                A(i, 2) = Psi3(particles, i + size * particles_set);
            }
        }
        det = A.determinant();
    }
    else if(m_particles == 2)
    {
        if(der_order == 0)
        {
            det = Psi1(particles, particles_set);
        }
        else if(der_order == 1)
        {
            det = GradiPsi1(particles, particles_set).at(grad_comp);
        }
        else if(der_order == 2)
        {
            det = LapliPsi1(particles, particles_set);
        }
    }
    if(m_particles == 12)
    {
        Eigen::MatrixXd A(size, size);
        for(int i = 0; i < size; i++)
        {
            if (i == row_changed)
            {
                if (der_order == 0)
                {
                    A(i, 0) = Psi1(particles, i + size * particles_set);
                    A(i, 1) = Psi2(particles, i + size * particles_set);
                    A(i, 2) = Psi3(particles, i + size * particles_set);
                    A(i, 3) = Psi4(particles, i + size * particles_set);
                    A(i, 4) = Psi5(particles, i + size * particles_set);
                    A(i, 5) = Psi6(particles, i + size * particles_set);
                }
                else if (der_order == 1)
                {
                    A(i, 0) = GradiPsi1(particles, i + size * particles_set).at(grad_comp);
                    A(i, 1) = GradiPsi2(particles, i + size * particles_set).at(grad_comp);
                    A(i, 2) = GradiPsi3(particles, i + size * particles_set).at(grad_comp);
                    A(i, 3) = GradiPsi4(particles, i + size * particles_set).at(grad_comp);
                    A(i, 4) = GradiPsi5(particles, i + size * particles_set).at(grad_comp);
                    A(i, 5) = GradiPsi6(particles, i + size * particles_set).at(grad_comp);

                }
                else if (der_order == 2)
                {
                    A(i, 0) = LapliPsi1(particles, i + size * particles_set);
                    A(i, 1) = LapliPsi2(particles, i + size * particles_set);
                    A(i, 2) = LapliPsi3(particles, i + size * particles_set);
                    A(i, 3) = LapliPsi4(particles, i + size * particles_set);
                    A(i, 4) = LapliPsi5(particles, i + size * particles_set);
                    A(i, 5) = LapliPsi6(particles, i + size * particles_set);
                }
            }
            else
            {
                A(i, 0) = Psi1(particles, i + size * particles_set);
                A(i, 1) = Psi2(particles, i + size * particles_set);
                A(i, 2) = Psi3(particles, i + size * particles_set);
                A(i, 3) = Psi4(particles, i + size * particles_set);
                A(i, 4) = Psi5(particles, i + size * particles_set);
                A(i, 5) = Psi6(particles, i + size * particles_set);
            }
        }
        det = A.determinant();
    }

    return det;
}


// Computes the sum of Laplacians of Psi_T divided by Psi_T for all particles in the system using the Fermion ansatz.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
double FermionNI::LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles)
{

    double det_up = SD(particles, 0, 0, 0, 0);      // Spin-up Slater Determinant.
    double det_down = SD(particles, 1, 0, 0, 0);    // Spin-down Slater Determinant.

    double sum = 0.0;
    for(int i = 0; i < m_particles / 2; i++)
    {
        // SD(particles, 0, i, 2, 0) is the Laplacian of the spin-up Slater determinant with respect to the i-th particle.
        sum += SD(particles, 0, i, 2, 0) / det_up;

        // SD(particles, 1, i, 2, 0) is the Laplacian of the spin-down Slater determinant with respect to the i-th particle.
        sum += SD(particles, 1, i, 2, 0) / det_down;
    }

    return sum;
}


// Evaluates the wave function for the given particle configuration using the fermionic ansatz.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the value of the wave function for the given particle configuration using the fermionic ansatz.
double FermionNI::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    // Spind-up Slater determinant times spin-down Slater determinant.
    return SD(particles, 0, 0, 0, 0) * SD(particles, 1, 0, 0, 0);
}


// Computes the sum of second derivatives of the wave function
// for the fermionic ansatz, which is required for calculating the local energy.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
double FermionNI::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    return LaplPsiTOverPsiT(particles);
}


// Computes the quantum force for a specific particle in the system using the fermionic ansatz.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int n - index of the particle for which to compute the quantum force.
//
// Output:  std::vector<double> - the quantum force for the specified particle.
std::vector<double> FermionNI::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    double det_up = SD(particles, 0, 0, 0, 0);      // Spin-up Slater Determinant.
    double det_down = SD(particles, 1, 0, 0, 0);    // Spin-down Slater Determinant.

    int half = m_particles / 2;  // Half the number of particles to differentiate spin-up and spin-down particles.
    std::vector<double> F(2);    // Initialize a vector to store the quantum force components.

    if (n < half)
    {
        // Compute the gradient of the spin-up Slater determinant with respect to the x-component and y-component.
        double grad_x = SD(particles, 0, n, 1, 0);
        double grad_y = SD(particles, 0, n, 1, 1);

        F[0] = 2.0 * grad_x / det_up;
        F[1] = 2.0 * grad_y / det_up;
    }
    else
    {
        int row = n - half;  // Adjust index for spin-down particles.

        // Compute the gradient of the spin-down Slater determinant with respect to the x-component and y-component.
        double grad_x = SD(particles, 1, row, 1, 0);
        double grad_y = SD(particles, 1, row, 1, 1);

        F[0] = 2.0 * grad_x / det_down;
        F[1] = 2.0 * grad_y / det_down;
    }

    return F;
}