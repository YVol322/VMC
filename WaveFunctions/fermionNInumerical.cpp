#include "fermionNInumerical.h"         // Include "fermionNInumerical" header file with declarations.

// Constructor of the FermionNInumerical class. Initializes the parameters for the fermion system.
// It sets the number of parameters, stores the variational parameter alpha, 
// and initializes the number of particles and angular frequency (omega).
//
// Input:   double alpha - variational parameter alpha;
//          int n_particles - number of particles in the system;
//          double omega - angular frequency of the system.
FermionNInumerical::FermionNInumerical(double alpha, int n_particles, double omega)
{
    assert(alpha >= 0);  // Ensure the variational parameter alpha is non-negative.

    m_numberOfParameters = 1;  // Only one parameter, alpha.

    m_parameters.reserve(m_numberOfParameters);  // Reserve space for the parameter.
    m_parameters.push_back(alpha);  // Store alpha as the variational parameter.

    m_particles = n_particles;  // Set the number of particles.

    m_omega = omega;  // Set the angular frequency.
    m_sqrt_om = sqrt(omega);  // Calculate the square root of omega.
}


// Fills the state vector x with particle positions from the system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  VectorXvar - the state vector containing the particle positions.
VectorXvar FermionNInumerical::fill_x(std::vector<std::unique_ptr<class Particle>>& particles)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    VectorXvar x(n_particles * n_dimensions);
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


// Computes the value of Psi1 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute Psi1.
//
// Output:  var - the value of Psi1 for the given particle.
var FermionNInumerical::psi1i(VectorXvar& x, int idx)
{
    var alpha = m_parameters.back();

    var x_ = x(idx);
    var y_ = x(idx + 1);
    var sum = x_ * x_ + y_ * y_;

    return exp(-alpha * m_omega * sum);
}


// Computes the value of Psi2 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute Psi2.
//
// Output:  var - the value of Psi2 for the given particle.
var FermionNInumerical::psi2i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx);
}


// Computes the value of Psi3 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute Psi3.
//
// Output:  var - the value of Psi3 for the given particle.
var FermionNInumerical::psi3i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx + 1);
}


// Computes the value of Psi4 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute Psi4.
//
// Output:  var - the value of Psi4 for the given particle.
var FermionNInumerical::psi4i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    return psi1 * x(idx) * x(idx + 1);
}


// Computes the value of Psi5 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute Psi5.
//
// Output:  var - the value of Psi5 for the given particle.
var FermionNInumerical::psi5i(VectorXvar& x, int idx)
{
    var x_ = x(idx);

    return psi1i(x, idx) * (x_ * x_ - 1);
}


// Computes the value of Psi6 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute Psi6.
//
// Output:  var - the value of Psi6 for the given particle.
var FermionNInumerical::psi6i(VectorXvar& x, int idx)
{
    var y_ = x(idx + 1);

    return psi1i(x, idx) * (y_ * y_ - 1);
}


// Computes the gradient of Psi1 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi1.
//
// Output:  VectorXvar - the gradient of Psi1 for the given particle.
VectorXvar FermionNInumerical::GradPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    // Using automatic differentiation to compute d psi1 / d x  and d psi1 / d y.
    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi1_x;
    grad(1) = psi1_y;

    return grad;
}


// Computes the gradient of Psi2 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi2.
//
// Output:  VectorXvar - the gradient of Psi2 for the given particle.
VectorXvar FermionNInumerical::GradPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    // Using automatic differentiation to compute d psi2 / d x  and d psi2 / d y.
    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi2_x;
    grad(1) = psi2_y;

    return grad;
}


// Computes the gradient of Psi3 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi3.
//
// Output:  VectorXvar - the gradient of Psi3 for the given particle.
VectorXvar FermionNInumerical::GradPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    // Using automatic differentiation to compute d psi3 / d x  and d psi3 / d y.
    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi3_x;
    grad(1) = psi3_y;

    return grad;
}


// Computes the gradient of Psi4 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi4.
//
// Output:  VectorXvar - the gradient of Psi4 for the given particle.
VectorXvar FermionNInumerical::GradPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    // Using automatic differentiation to compute d psi4 / d x  and d psi4 / d y.
    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi4_x;
    grad(1) = psi4_y;

    return grad;
}


// Computes the gradient of Psi5 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi5.
//
// Output:  VectorXvar - the gradient of Psi5 for the given particle.
VectorXvar FermionNInumerical::GradPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    // Using automatic differentiation to compute d psi5 / d x  and d psi5 / d y.
    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi5_x;
    grad(1) = psi5_y;

    return grad;
}


// Computes the gradient of Psi6 for the specified particle using the state vector x.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi6.
//
// Output:  VectorXvar - the gradient of Psi6 for the given particle.
VectorXvar FermionNInumerical::GradPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    // Using automatic differentiation to compute d psi6 / d x  and d psi6 / d y.
    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi6_x;
    grad(1) = psi6_y;

    return grad;
}


// Computes the Laplacian of Psi1 for the specified particle using the state vector `x` and automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi1.
//
// Output:  var - the Laplacian of Psi1 for the given particle.
var FermionNInumerical::LaplPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    // Using automatic differentiation to compute d^2 psi1 / d x^2  and d^2 psi1 / d y^2.
    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    auto [psi1_xx] = derivativesx(psi1_x, wrt(x(idx)));
    auto [psi1_yy] = derivativesx(psi1_y, wrt(x(idx + 1)));

    return psi1_xx + psi1_yy;
}


// Computes the Laplacian of Psi2 for the specified particle using the state vector `x` and automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi2.
//
// Output:  var - the Laplacian of Psi2 for the given particle.
var FermionNInumerical::LaplPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    // Using automatic differentiation to compute d^2 psi2 / d x^2  and d^2 psi2 / d y^2.
    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    auto [psi2_xx] = derivativesx(psi2_x, wrt(x(idx)));
    auto [psi2_yy] = derivativesx(psi2_y, wrt(x(idx + 1)));

    return psi2_xx + psi2_yy;
}


// Computes the Laplacian of Psi3 for the specified particle using the state vector `x` and automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi3.
//
// Output:  var - the Laplacian of Psi3 for the given particle.
var FermionNInumerical::LaplPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    // Using automatic differentiation to compute d^2 psi3 / d x^2  and d^2 psi3 / d y^2.
    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    auto [psi3_xx] = derivativesx(psi3_x, wrt(x(idx)));
    auto [psi3_yy] = derivativesx(psi3_y, wrt(x(idx + 1)));

    return psi3_xx + psi3_yy;
}


// Computes the Laplacian of Psi4 for the specified particle using the state vector `x` and automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi4.
//
// Output:  var - the Laplacian of Psi4 for the given particle.
var FermionNInumerical::LaplPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    // Using automatic differentiation to compute d^2 psi4 / d x^2  and d^2 psi4 / d y^2.
    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    auto [psi4_xx] = derivativesx(psi4_x, wrt(x(idx)));
    auto [psi4_yy] = derivativesx(psi4_y, wrt(x(idx + 1)));

    return psi4_xx + psi4_yy;
}


// Computes the Laplacian of Psi5 for the specified particle using the state vector `x` and automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi5.
//
// Output:  var - the Laplacian of Psi5 for the given particle.
var FermionNInumerical::LaplPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    // Using automatic differentiation to compute d^2 psi5 / d x^2  and d^2 psi5 / d y^2.
    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    auto [psi5_xx] = derivativesx(psi5_x, wrt(x(idx)));
    auto [psi5_yy] = derivativesx(psi5_y, wrt(x(idx + 1)));

    return psi5_xx + psi5_yy;
}


// Computes the Laplacian of Psi6 for the specified particle using the state vector `x` and automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi6.
//
// Output:  var - the Laplacian of Psi6 for the given particle.
var FermionNInumerical::LaplPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    // Using automatic differentiation to compute d^2 psi6 / d x^2  and d^2 psi6 / d y^2.
    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    auto [psi6_xx] = derivativesx(psi6_x, wrt(x(idx)));
    auto [psi6_yy] = derivativesx(psi6_y, wrt(x(idx + 1)));

    return psi6_xx + psi6_yy;
}


// Computes the Slater determinant or its derivatives (gradient or Laplacian) for the fermionic system using the given parameters.
//
// Input:   VectorXvar& x - state vector containing particle positions;
//          int particles - number of particles in the system;
//          int row_changed - row index of the particle that has been changed;
//          int der_order - order of the derivatives to compute (0 - value, 1 - gradient, 2 - Laplacian);
//          int grad_comp - gradient component to compute (if applicable, typically x or y component).
//
// Output:  double - the computed Slater determinant value or its derivative (gradient or Laplacian).
double FermionNInumerical::SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp)
{
    int size = m_particles / 2;

    Eigen::MatrixXd A(size, size);
    double det = 0;
    if(m_particles == 6)
    {
        for(int i = 0; i < size; i++)
        {
            int idx = i * 2 + m_particles * particles;
            {
                if (i == row_changed)
                {
                    if (der_order == 0)
                    {
                        A(i, 0) = val(psi1i(x, idx));
                        A(i, 1) = val(psi2i(x, idx));
                        A(i, 2) = val(psi3i(x, idx));
                    }
                    else if (der_order == 1)
                    {
                        A(i, 0) = val((GradPsi1i(x, idx))(grad_comp));
                        A(i, 1) = val((GradPsi2i(x, idx))(grad_comp));
                        A(i, 2) = val((GradPsi3i(x, idx))(grad_comp));
                    }
                    else if (der_order == 2)
                    {
                        A(i, 0) = val(LaplPsi1i(x, idx));
                        A(i, 1) = val(LaplPsi2i(x, idx));
                        A(i, 2) = val(LaplPsi3i(x, idx));
                    }
                }
                else
                {
                    A(i, 0) = val(psi1i(x, idx));
                    A(i, 1) = val(psi2i(x, idx));
                    A(i, 2) = val(psi3i(x, idx));
                }
            }
        }
        det = A.determinant();
    }
    else if(m_particles == 2)
    {
        int idx = m_particles * particles;
        if(der_order == 0)
        {
            det = val(psi1i(x, idx));
        }
        else if(der_order == 1)
        {
            det = val((GradPsi1i(x, idx))(grad_comp));
        }
        else if(der_order == 2)
        {
            det = val(LaplPsi1i(x, idx));
        }
    }
    else if(m_particles == 12)
    {
        for(int i = 0; i < size; i++)
        {
            int idx = i * 2 + m_particles * particles;
            {
                if (i == row_changed)
                {
                    if (der_order == 0)
                    {
                        A(i, 0) = val(psi1i(x, idx));
                        A(i, 1) = val(psi2i(x, idx));
                        A(i, 2) = val(psi3i(x, idx));
                        A(i, 3) = val(psi4i(x, idx));
                        A(i, 4) = val(psi5i(x, idx));
                        A(i, 5) = val(psi6i(x, idx));
                    }
                    else if (der_order == 1)
                    {
                        A(i, 0) = val((GradPsi1i(x, idx))(grad_comp));
                        A(i, 1) = val((GradPsi2i(x, idx))(grad_comp));
                        A(i, 2) = val((GradPsi3i(x, idx))(grad_comp));
                        A(i, 3) = val((GradPsi4i(x, idx))(grad_comp));
                        A(i, 4) = val((GradPsi5i(x, idx))(grad_comp));
                        A(i, 5) = val((GradPsi6i(x, idx))(grad_comp));
                    }
                    else if (der_order == 2)
                    {
                        A(i, 0) = val(LaplPsi1i(x, idx));
                        A(i, 1) = val(LaplPsi2i(x, idx));
                        A(i, 2) = val(LaplPsi3i(x, idx));
                        A(i, 3) = val(LaplPsi4i(x, idx));
                        A(i, 4) = val(LaplPsi5i(x, idx));
                        A(i, 5) = val(LaplPsi6i(x, idx));
                    }
                }
                else
                {
                    A(i, 0) = val(psi1i(x, idx));
                    A(i, 1) = val(psi2i(x, idx));
                    A(i, 2) = val(psi3i(x, idx));
                    A(i, 3) = val(psi4i(x, idx));
                    A(i, 4) = val(psi5i(x, idx));
                    A(i, 5) = val(psi6i(x, idx));
                }
            }
        }
        det = A.determinant();
    }

    return det;
}


// Computes the sum of Laplacians of Psi_T divided by Psi_T for all particles in the system using the Fermion ansatz.
//
// Input:   VectorXvar& x - state vector containing particle positions.
//
// Output:  double - the computed sum of Laplacians of Psi_T divided by Psi_T for all particles in the system.
double FermionNInumerical::LaplPsiTOverPsiT(VectorXvar& x)
{
    double Psi_up = SD(x, 0, 0, 0, 0);         // Spin-up Slater Determinant.
    double Psi_down = SD(x, 1, 0, 0, 0);       // Spin-down Slater Determinant.

    double sum = 0.0;
    for (int i = 0; i < m_particles; i++)
    {
        if (i < m_particles / 2)
        {
            // SD(particles, 0, i, 2, 0) is the Laplacian of the spin-up Slater determinant with respect to the i-th particle.
            sum += SD(x, 0, i, 2, 0) * Psi_down;
        }
        else
        {
            int j = i - m_particles / 2;

            // SD(particles, 1, j, 2, 0) is the Laplacian of the spin-down Slater determinant with respect to the i-th particle.
            sum += Psi_up * SD(x, 1, j, 2, 0);
        }
    }

    return sum / (Psi_up * Psi_down);
}


// Computes the wave function for the fermionic system by evaluating the Slater determinants for spin-up and spin-down particles.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed wave function (Psi) for the given particle configuration.
double FermionNInumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    double Psi_up = SD(x, 0, 0, 0, 0);         // Spin-up Slater Determinant.
    double Psi_down = SD(x, 1, 0, 0, 0);       // Spin-down Slater Determinant.

    return Psi_up * Psi_down;
}


// Computes the sum of second derivatives (Laplacian operator) for the fermionic system, which is needed to compute the local energy.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of second derivatives for the given particle configuration.
double FermionNInumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    return LaplPsiTOverPsiT(x);
}


// Computes the quantum force for a specific particle in the fermionic system using the fermionic ansatz.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int n - index of the particle for which to compute the quantum force.
//
// Output:  std::vector<double> - the quantum force for the specified particle.
std::vector<double> FermionNInumerical::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    VectorXvar x = fill_x(particles);
    int size = m_particles / 2;             // Half the number of particles to differentiate spin-up and spin-down particles.


    double det_up = SD(x, 0, 0, 0, 0);      // Spin-up Slater Determinant.
    double det_down = SD(x, 1, 0, 0, 0);    // Spin-down Slater Determinant.

    std::vector<double> F;                  // Initialize a vector to store the quantum force components.

    if(n < m_particles/2)
    {
        // Compute the gradient of the spin-up Slater determinant with respect to the x-component and y-component.
        double grad_x = SD(x, 0, n, 1, 0);
        double grad_y = SD(x, 0, n, 1, 1);

        F = {2 * grad_x / det_up, 2 * grad_y / det_up};
    }
    else
    {
        // Compute the gradient of the spin-down Slater determinant with respect to the x-component and y-component.
        double grad_x = SD(x, 1, n - size, 1, 0);
        double grad_y = SD(x, 1, n - size, 1, 1);

        F = {2 * grad_x / det_down, 2 * grad_y / det_down};
    }

    return F;
}
