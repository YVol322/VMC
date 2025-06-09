#include "fermion.h"        // Include "fermion" header file with declarations.


// Constructor of the Fermion class. Initializes the parameters for the fermion system.
//
// Input:   double alpha - variational parameter alpha;
//          std::vector<double> beta - vector of variational parameters for Jastrow or Pade-Jastrow;
//          int mode - 0 for Jastrow ansatz, 1 for Pade-Jastrow ansatz;
//          int n_particles - number of particles in the system;
//          double omega - angular frequency of the system.
Fermion::Fermion(double alpha, std::vector<double> beta, int mode, int n_particles, double omega)
{

    int n_betas = beta.size();
    m_numberOfParameters = n_betas + 1;         // n_betas beta parameters and one alpha parameter.

    m_parameters.reserve(m_numberOfParameters);     // Reserve space for the parameter.
    m_parameters.insert(m_parameters.end(), beta.begin(), beta.end());  // First store beta as the variational parameters.
    m_parameters.push_back(alpha);  // Then store alpha as the variational parameter.

    m_particles = n_particles;      // Set the number of particles.
    m_mode = mode;                  // Set the mode: 0 - Jastrow ansatz, 1 - Pade-Jastrow ansatz.
    m_omega = omega;                // Set the angular frequency.
    m_sqrt_om = sqrt(omega);        // Calculate the square root of omega.
}


// Computes the Jastrow factor for the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed Jastrow factor for the given particle configuration.
double Fermion::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;


    // These loops go as (i,j) = (0,1), (0,2), ..., (0, n), (1, 2), ..., (1, n), ..., (n-1, n).
    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {
            double rij = m_sqrt_om * r_ij(particles, i, j);

            // Transform 2D index (i, j) into 1D index.
            int idx = BetaIndex(i,j);
            double beta_ij = m_parameters[idx];

            sum += beta_ij * rij;
        }
    }

    return exp(sum);
}


// Computes the Pade-Jastrow factor for the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed Pade-Jastrow factor for the given particle configuration.
double Fermion::PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;
    double beta = m_parameters[0];

    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {
            // Using the convention that particles with indices 0, 1, ..., n_particles / 2 are spin-up, and the rest are spin-down, 
            // we compute the constant 'a' for the given pair of particles. Since j is always greater than i, we only need to check 
            // if i is spin-up (i < n_particles / 2) and j is spin-down (j >= n_particles / 2) to compute 'a'.
            double aij = a_ij(i, j);
            double rij = m_sqrt_om * r_ij(particles, i, j);

            sum += aij * rij / (1 + beta * rij);
        }
    }

    return exp(sum);
}


// Computes the value of Psi1 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi1.
//
// Output:  double - the value of Psi1 for the given particle.
double Fermion::Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    double r2 = r_squared(particles, part_idx);

    return exp(-alpha * r2 * m_omega);
}


// Computes the value of Psi2 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute Psi2.
//
// Output:  double - the value of Psi2 for the given particle.
double Fermion::Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
double Fermion::Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
double Fermion::Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
double Fermion::Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
double Fermion::Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
std::vector<double> Fermion::GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
std::vector<double> Fermion::GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
std::vector<double> Fermion::GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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
std::vector<double> Fermion::GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi4 = {y * psi1 * (1 - 2 * alpha * m_omega * x * x),
                                        x * psi1 * (1 - 2 * alpha * m_omega *y * y)};

    return grad_psi4;
}


// Computes the gradient of Psi5 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi5.
//
// Output:  std::vector<double> - the gradient of Psi5 for the given particle.
std::vector<double> Fermion::GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi5 = {
        2.0 * x * psi1 * (1.0 - alpha * m_omega * x * x + alpha * m_omega),
       -2.0 * alpha * m_omega * y * (x * x - 1.0) * psi1
    };


    return grad_psi5;
}


// Computes the gradient of Psi6 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi6.
//
// Output:  std::vector<double> - the gradient of Psi6 for the given particle.
std::vector<double> Fermion::GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi6 = {
        -2.0 * alpha * m_omega * x * (y * y - 1.0) * psi1,
         2.0 * y * psi1 * (1.0 - alpha * m_omega * y * y + alpha * m_omega)
    };

    return grad_psi6;
}


// Computes the Laplacian of Psi1 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi1.
//
// Output:  double - the Laplacian of Psi1 for the given particle.
double Fermion::LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-4 * alpha * m_omega + 4 * alpha * alpha * m_omega * m_omega * (x * x + y * y) ) * Psi1(particles, part_idx);
}


// Computes the Laplacian of Psi2 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi2.
//
// Output:  double - the Laplacian of Psi2 for the given particle.
double Fermion::LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-8 * alpha * m_omega * x + 4 * alpha * alpha * m_omega * m_omega * 
                                                                x * (x * x + y * y) ) * Psi1(particles, part_idx);
}


// Computes the Laplacian of Psi3 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi3.
//
// Output:  double - the Laplacian of Psi3 for the given particle.
double Fermion::LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-8 * alpha * m_omega * y + 4 * alpha * alpha * m_omega * m_omega * 
                                                                y * (x * x + y * y) ) * Psi1(particles, part_idx);
}


// Computes the Laplacian of Psi4 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi4.
//
// Output:  double - the Laplacian of Psi4 for the given particle.
double Fermion::LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return x * y * Psi1(particles, part_idx) * (4 * alpha * alpha * m_omega * m_omega * (x * x + y * y)
                                                                                        - 12 * alpha * m_omega);
}


// Computes the Laplacian of Psi5 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi5.
//
// Output:  double - the Laplacian of Psi5 for the given particle.
double Fermion::LapliPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    double x2 = x * x;
    double y2 = y * y;
    double common = x2 - 1;

    double laplacian_psi5 = 2 * psi1 *
    (-4 * alpha * m_omega * x2 +
    alpha * m_omega * common * (2 * alpha * m_omega * x2 - 1)+
    alpha * m_omega * common * (2 * alpha * m_omega * y2 - 1) + 1);

    return laplacian_psi5;
}


// Computes the Laplacian of Psi6 for the specified particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of Psi6.
//
// Output:  double - the Laplacian of Psi6 for the given particle.
double Fermion::LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    double x2 = x * x;
    double y2 = y * y;
    double common = y2 - 1;

    double laplacian_psi6 = 2.0 * psi1 * (
        -4.0 * alpha * m_omega * y2
      + (alpha * m_omega) * common * (2.0 * alpha * m_omega * x2 - 1.0)
      + (alpha * m_omega) * common * (2.0 * alpha * m_omega * y2 - 1.0)
      + 1.0
    );

    return laplacian_psi6;
}


// Computes the Slater determinant or its gradient component or Laplacian depending on the derivative order.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int particles_set - set of particles: 0 for spin-up, 1 for spin-down;
//          int row_changed - row index of the particle that was changed;
//          int der_order - order of the derivatives to compute: 0 for Psi, 1 for gradient, 2 for Laplacian;
//          int grad_comp - component of the gradient to compute (if applicable).
//
// Output:  double - the computed value based on the requested derivative order (Slater determinant, gradient, or Laplacian).
double Fermion::SD
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


// Computes the gradient of the Jastrow factor for a given particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of the Jastrow factor.
//
// Output:  std::vector<double> - the gradient of the Jastrow factor for the given particle.
std::vector<double> Fermion::GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    Particle& particle_i = *(particles[part_idx]);

    std::vector<double> grad(n_dimensions, 0.0);

    for(int j = 0; j < n_particles; j++)
    {
        // Skip the interaction computation for the particle with itself (i.e., if j == part_idx, the particle doesn't interact with itself).
        if(j == part_idx) continue;

        Particle& particle_j = *(particles[j]);
        double ril = r_ij(particles, part_idx, j);

        // Extracting the correct beta parameter index for the interaction between particle 'part_idx' and particle 'j'.
        int idx = BetaIndex(part_idx, j);
        double beta_il = m_parameters[idx];

        for(int d = 0; d < n_dimensions; d++)
        {
            double coord_i = particle_i.getPosition()[d];
            double coord_j = particle_j.getPosition()[d];

            grad[d] += beta_il * m_sqrt_om *(coord_i - coord_j) / ril;
        }
    }

    return grad;
}


// Computes the gradient of the Pade-Jastrow factor for a given particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_inx - index of the particle for which to compute the gradient of the Pade-Jastrow factor.
//
// Output:  std::vector<double> - the gradient of the Pade-Jastrow factor for the given particle.
std::vector<double> Fermion::GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    int n_particles = particles.size();
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    double beta = m_parameters.at(0);

    Particle particle_i = *(particles[part_idx]);

    std::vector<double> grad(n_dimensions, 0.0);

    for(int j = 0; j < n_particles; j++)
    {
        // Skip the interaction computation for the particle with itself (i.e., if j == part_idx, the particle doesn't interact with itself).
        if(j == part_idx) continue;

        Particle particle_j = *(particles[j]);
        double ril = r_ij(particles, part_idx, j);

        // Computing the constant 'a' for the interaction between particles 'part_idx' and 'j'.
        double aij = a_ij(part_idx, j);
        double denom = (1 + beta * m_sqrt_om *ril);

        for(int d = 0; d < n_dimensions; d++)
        {
            double coord_i = particle_i.getPosition()[d];
            double coord_j = particle_j.getPosition()[d];

            grad[d] += aij * m_sqrt_om * (coord_i - coord_j) / (ril * denom * denom);
        }
    }

    return grad;
}


// Computes the Laplacian of the Jastrow factor for a given particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the Laplacian of the Jastrow factor.
//
// Output:  double - the Laplacian of the Jastrow factor for the given particle.
double Fermion::LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    std::vector<double> grad_i = GradiJOverJ(particles, part_idx);

    // Compute the squared gradient of the logarithm of the Jastrow factor (Grad^2 ln(J)) term.
    // This sums the squared components of the gradient vector grad_i for each dimension.
    double sum1 = 0.0;
    for (int d = 0; d < n_dimensions; d++)
    {
        sum1 += grad_i[d] * grad_i[d];
    }

    // Compute the Laplacian of the logarithm of the Jastrow factor (Lapl(ln(J))) term.
    // This calculates the contribution from each particle j (except the current particle 'part_idx')
    // based on the relative distance between particles 'part_idx' and 'j', weighted by the corresponding beta parameter.
    double sum2 = 0.0;
    for (int j = 0; j < n_particles; j++)
    {
        // Skip the interaction computation for the particle with itself (i.e., if j == part_idx, the particle doesn't interact with itself).
        if (j == part_idx) continue;

        // Extracting the correct beta parameter index for the interaction between particle 'part_idx' and particle 'j'.
        int idx = BetaIndex(part_idx, j);
        double beta_il = m_parameters[idx];

        double ril = r_ij(particles, part_idx, j);

        sum2 += m_sqrt_om * beta_il / ril;
    }

    return sum1 + sum2;
}


// Computes the Laplacian of the Pade-Jastrow factor for a given particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_inx - index of the particle for which to compute the Laplacian of the Pade-Jastrow factor.
//
// Output:  double - the Laplacian of the Pade-Jastrow factor for the given particle.
double Fermion::LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    std::vector<double> grad_i = GradiPJOverPJ(particles, part_idx);

    // Compute the squared gradient of the logarithm of the Pade-Jastrow factor (Grad^2 ln(P)) term.
    // This sums the squared components of the gradient vector grad_i for each dimension.
    double sum1 = 0.0;
    for (int d = 0; d < n_dimensions; d++)
    {
        sum1 += grad_i[d] * grad_i[d];
    }

    // Compute the Laplacian of the logarithm of the Pade-Jastrow factor (Lapl(ln(P))) term.
    // This calculates the contribution from each particle j (except the current particle 'part_idx')
    // based on the relative distance between particles 'part_idx' and 'j', weighted by the corresponding beta parameter.
    double sum2 = 0.0;

    // Get the single variational parameter beta, which is stored in the first position of the parameters vector.
    double beta = m_parameters[0];
    for (int j = 0; j < n_particles; j++)
    {
        if (j == part_idx) continue;

        // Computing the constant 'a' for the interaction between particles 'part_idx' and 'j'.
        double aij = a_ij(part_idx, j);
        double ril = r_ij(particles, part_idx, j);
        double t = 1 + m_sqrt_om * beta * ril;

        sum2 += aij * m_sqrt_om * (1/(ril * t * t) - 2 * m_sqrt_om * beta / (t * t * t));
    }

    return sum1 + sum2;
}


// Computes the sum of Laplacians of Psi_T divided by Psi_T for all particles in the system using the Fermion ansatz.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
double Fermion::LaplPsiTOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double det_up = SD(particles, 0, 0, 0, 0);         // Spin-up Slater Determinant.
    double det_down = SD(particles, 1, 0, 0, 0);       // Spin-down Slater Determinant.

    double sum = 0.0;
    for(int i = 0; i < m_particles / 2; i++)
    {
        // SD(particles, 0, i, 2, 0) is the Laplacian of the spin-up Slater determinant with respect to the i-th particle.
        sum += SD(particles, 0, i, 2, 0) / det_up;

        // SD(particles, 1, j, 2, 0) is the Laplacian of the spin-down Slater determinant with respect to the i-th particle.
        sum += SD(particles, 1, i, 2, 0) / det_down;
    }

    return sum;
}


// Computes the gradient of Psi_T multiplied by the gradient of the Jastrow factor (or the Pade-Jastrow factor) 
// over Psi_T for a given particle.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          double part_idx - index of the particle for which to compute the gradient of Psi_T multiplied by
//                            gradient of J (or PJ) over Psi_T.
//
// Output:  double - the computed value of the gradient of Psi_T multiplied by gradient of J (or PJ) over Psi_T.
double Fermion::GradiPsiTGradiJOverPsiT(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    // If m_mode is 0, compute the gradient of J (Jastrow factor), otherwise compute the gradient of Pade-Jastrow (PJ) factor.
    std::vector<double> grad_J = (m_mode == 0) ? GradiJOverJ(particles, part_idx) : GradiPJOverPJ(particles, part_idx);

    double det_up = SD(particles, 0, 0, 0, 0);         // Spin-up Slater Determinant.
    double det_down = SD(particles, 1, 0, 0, 0);       // Spin-down Slater Determinant.

    // If part_idx is less than half of m_particles, compute only its contribution to the gradient, 
    // since the other Slater determinant cancels out due to symmetry.
    double sum = 0.0;
    if (part_idx < m_particles/2)
    {
        for (int i = 0; i < n_dimensions; i++)
        {
            // SD(particles, 0, part_idx, 1, i) computes the i-th component of the gradient of the spin-up Slater determinant 
            // with respect to the part_idx-th particle. The first argument '0' refers to spin-up particles, 
            // 'part_idx' specifies the particle for which the gradient is calculated, '1' indicates that we are computing the gradient (first derivative),
            // and 'i' specifies the component of the gradient (e.g., x or y direction).
            sum += SD(particles, 0, part_idx, 1, i) / det_up * grad_J[i];
        }
    }
    else
    {
        int idx = part_idx - m_particles/2; // Adjusts the part_idx to access elements of the spin-down determinant correctly. 
        for (int i = 0; i < n_dimensions; i++)
        {
            // SD(particles, 1, part_idx, 1, i) computes the i-th component of the gradient of the spin-down Slater determinant 
            // with respect to the part_idx-th particle. The first argument '1' refers to spin-down particles, 
            // 'part_idx' specifies the particle for which the gradient is calculated, '1' indicates that we are computing the gradient (first derivative),
            // and 'i' specifies the component of the gradient (e.g., x or y direction).
            sum += SD(particles, 1, idx, 1, i) / det_down * grad_J[i];
        }
    }

    return 2.0 * sum;
}


// WaveFunction class function that evaluates the wave function for the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed value of the wave function for the given particle configuration.
double Fermion::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    // If m_mode is 0, compute the gradient of J (Jastrow factor), otherwise compute the gradient of Pade-Jastrow (PJ) factor.
    double J = (m_mode == 0) ? Jastrow(particles) : PadeJastrow(particles);

    return SD(particles, 0, 0, 0, 0) * SD(particles, 1, 0, 0, 0) * J;
}


// WaveFunction class function that computes the sum of second derivatives (Laplace operator) of the wave function
// for the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
double Fermion::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double LaplJ = 0;
    double GradGrad = 0;

    // This block of code computes the sum of Laplacians and the sum of gradient terms for each particle
    // depending on the Jastrow ansatz mode. 
    if(m_mode == 0)
    {
        for(int i = 0; i < m_particles; i++)
        {
            LaplJ += LapliJOverJ(particles, i);
            GradGrad += GradiPsiTGradiJOverPsiT(particles, i);
        }
    }
    else
    {
        for(int i = 0; i < m_particles; i++)
        {
            LaplJ += LapliPJOverPJ(particles, i);
            GradGrad += GradiPsiTGradiJOverPsiT(particles, i);
        }
    }

    // Add the sum of Laplacians of the ansatz (Psi_T) divided by Psi_T .
    return LaplPsiTOverPsiT(particles) + LaplJ + GradGrad;
}


// WaveFunction class function that computes the quantum force for a specific particle in the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int i - index of the particle for which to compute the quantum force.
//
// Output:  std::vector<double> - the quantum force for the specified particle.
std::vector<double> Fermion::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int i)
{
    int half = m_particles/2;
    int spinBlock = (i < half ? 0 : 1);         // Set the spin block: 0 for spin-up, 1 for spin-down.
    int localIdx  = (i < half ? i : i - half);  // Adjust the index for spin-down particles.

    // Compute the gradient of the Slater determinant for the x and y components (dlnD_dx and dlnD_dy).
    // SD(particles, spinBlock, localIdx, 1, 0) computes the gradient w.r.t. x, and SD(particles, spinBlock, 0, 0, 0) is the determinant.
    double dlnD_dx = SD(particles, spinBlock, localIdx, 1, 0) / SD(particles, spinBlock, 0, 0, 0);
    double dlnD_dy = SD(particles, spinBlock, localIdx, 1, 1) / SD(particles, spinBlock, 0, 0, 0);

    // Compute the gradient of the Jastrow or Pade-Jastrow factor based on the mode selected (m_mode).
    auto dlnJ = (m_mode == 0) ? GradiJOverJ(particles, i) : GradiPJOverPJ(particles, i);

    std::vector<double> F = { 2*(dlnD_dx + dlnJ[0]), 2*(dlnD_dy + dlnJ[1]) };
    
    return F;
}