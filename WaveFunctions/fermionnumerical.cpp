#include "fermionnumerical.h"       // Include "fermionnumerical" header file with declarations.


// Constructor of the FermionNumerical class. Initializes the parameters for the fermion system.
// Input:   double alpha - variational parameter alpha;
//          std::vector<double> beta - vector of variational parameters for Jastrow or Pade-Jastrow;
//          int mode - 0 for Jastrow ansatz, 1 for Pade-Jastrow ansatz;
//          int n_particles - number of particles in the system;
//          double omega - angular frequency of the system.
FermionNumerical::FermionNumerical(double alpha, std::vector<double>beta, int mode, int n_particles, double omega)
{
    assert(alpha >= 0);

    int n_betas = beta.size();
    m_numberOfParameters = n_betas + 1;

    m_parameters.reserve(m_numberOfParameters);
    m_parameters.insert(m_parameters.end(), beta.begin(), beta.end());
    m_parameters.push_back(alpha);

    m_particles = n_particles;
    m_mode = mode;

    m_omega = omega;
    m_sqrt_om = sqrt(omega);
}


// Fills the state vector x with the particle positions for automatic differentiation.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  VectorXvar - state vector containing the particle positions.
VectorXvar FermionNumerical::fill_x(std::vector<std::unique_ptr<class Particle>>& particles)
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


// Computes the Jastrow factor for the fermionic system using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions.
//
// Output:  var - the computed Jastrow factor for the given particle configuration.
var FermionNumerical::Jastrow(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;

    var sum1 = 0.0;
    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {
            var sum2 = 0.0;
            for (int d = 0; d < n_dimensions; d++)
            {
                var delta = x(i * n_dimensions + d) - x(j * n_dimensions + d);
                sum2 += delta * delta;
            }
            var rij = m_sqrt_om * sqrt(sum2);

            int idx = BetaIndex(i,j);
            var beta_ij = m_parameters[idx];

            sum1 += beta_ij * rij;
        }
    }

    return exp(sum1);
}


// Computes the Pade-Jastrow factor for the fermionic system using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions.
//
// Output:  var - the computed Pade-Jastrow factor for the given particle configuration.
var FermionNumerical::PadeJastrow(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;

    var beta = m_parameters[0];
    var sum1 = 0.0;
    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {

            var sum2 = 0.0;
            for (int d = 0; d < n_dimensions; d++)
            {
                var delta = x(i * n_dimensions + d) - x(j * n_dimensions + d);
                sum2 += delta * delta;
            }
            var rij = m_sqrt_om * sqrt(sum2);

            var aij = a_ij(i, j);

            sum1 += aij * rij / (1 + beta * rij);
        }
    }

    return exp(sum1);
}


// Computes the value of Psi1 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute Psi1.
//
// Output:  var - the value of Psi1 for the given particle.
var FermionNumerical::psi1i(VectorXvar& x, int idx)
{
    var alpha = m_parameters.back();

    var x_ = x(idx);
    var y_ = x(idx + 1);
    var sum = x_ * x_ + y_ * y_;

    return exp(-alpha * m_omega * sum);
}


// Computes the value of Psi2 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute Psi2.
//
// Output:  var - the value of Psi2 for the given particle.
var FermionNumerical::psi2i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx);
}


// Computes the value of Psi3 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute Psi3.
//
// Output:  var - the value of Psi3 for the given particle.
var FermionNumerical::psi3i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx + 1);
}


// Computes the value of Psi4 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute Psi4.
//
// Output:  var - the value of Psi4 for the given particle.
var FermionNumerical::psi4i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    return psi1 * x(idx) * x(idx + 1);
}


// Computes the value of Psi5 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute Psi5.
//
// Output:  var - the value of Psi5 for the given particle.
var FermionNumerical::psi5i(VectorXvar& x, int idx)
{
    var x_ = x(idx);

    return psi1i(x, idx) * (x_ * x_ - 1);
}


// Computes the value of Psi6 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute Psi6.
//
// Output:  var - the value of Psi6 for the given particle.
var FermionNumerical::psi6i(VectorXvar& x, int idx)
{
    var y_ = x(idx + 1);

    return psi1i(x, idx) * (y_ * y_ - 1);
}


// Computes the gradient of Psi1 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi1.
//
// Output:  VectorXvar - the gradient of Psi1 for the given particle.
VectorXvar FermionNumerical::GradPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    // Using automatic differentiation to compute d psi1 / d x  and d psi1 / d y.
    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi1_x;
    grad(1) = psi1_y;

    return grad;
}


// Computes the gradient of Psi2 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi2.
//
// Output:  VectorXvar - the gradient of Psi2 for the given particle.
VectorXvar FermionNumerical::GradPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    // Using automatic differentiation to compute d psi2 / d x  and d psi2 / d y.
    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi2_x;
    grad(1) = psi2_y;

    return grad;
}


// Computes the gradient of Psi3 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi3.
//
// Output:  VectorXvar - the gradient of Psi3 for the given particle.
VectorXvar FermionNumerical::GradPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    // Using automatic differentiation to compute d psi3 / d x  and d psi3 / d y.
    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi3_x;
    grad(1) = psi3_y;

    return grad;
}


// Computes the gradient of Psi4 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi4.
//
// Output:  VectorXvar - the gradient of Psi4 for the given particle.
VectorXvar FermionNumerical::GradPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    // Using automatic differentiation to compute d psi4 / d x  and d psi4 / d y.
    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi4_x;
    grad(1) = psi4_y;

    return grad;
}


// Computes the gradient of Psi5 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi5.
//
// Output:  VectorXvar - the gradient of Psi5 for the given particle.
VectorXvar FermionNumerical::GradPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    // Using automatic differentiation to compute d psi5 / d x  and d psi5 / d y.
    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi5_x;
    grad(1) = psi5_y;

    return grad;
}


// Computes the gradient of Psi6 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of Psi6.
//
// Output:  VectorXvar - the gradient of Psi6 for the given particle.
VectorXvar FermionNumerical::GradPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    // Using automatic differentiation to compute d psi6 / d x  and d psi6 / d y.
    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi6_x;
    grad(1) = psi6_y;

    return grad;
}


// Computes the Laplacian of Psi1 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi1.
//
// Output:  var - the Laplacian of Psi1 for the given particle.
var FermionNumerical::LaplPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    // Using automatic differentiation to compute d^2 psi1 / d x^2  and d^2 psi1 / d y^2.
    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    auto [psi1_xx] = derivativesx(psi1_x, wrt(x(idx)));
    auto [psi1_yy] = derivativesx(psi1_y, wrt(x(idx + 1)));

    return psi1_xx + psi1_yy;
}


// Computes the Laplacian of Psi2 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi2.
//
// Output:  var - the Laplacian of Psi2 for the given particle.
var FermionNumerical::LaplPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    // Using automatic differentiation to compute d^2 psi2 / d x^2  and d^2 psi2 / d y^2.
    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    auto [psi2_xx] = derivativesx(psi2_x, wrt(x(idx)));
    auto [psi2_yy] = derivativesx(psi2_y, wrt(x(idx + 1)));

    return psi2_xx + psi2_yy;
}


// Computes the Laplacian of Psi3 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi3.
//
// Output:  var - the Laplacian of Psi3 for the given particle.
var FermionNumerical::LaplPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    // Using automatic differentiation to compute d^2 psi3 / d x^2  and d^2 psi3 / d y^2.
    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    auto [psi3_xx] = derivativesx(psi3_x, wrt(x(idx)));
    auto [psi3_yy] = derivativesx(psi3_y, wrt(x(idx + 1)));

    return psi3_xx + psi3_yy;
}


// Computes the Laplacian of Psi4 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi4.
//
// Output:  var - the Laplacian of Psi4 for the given particle.
var FermionNumerical::LaplPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    // Using automatic differentiation to compute d^2 psi4 / d x^2  and d^2 psi4 / d y^2.
    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    auto [psi4_xx] = derivativesx(psi4_x, wrt(x(idx)));
    auto [psi4_yy] = derivativesx(psi4_y, wrt(x(idx + 1)));

    return psi4_xx + psi4_yy;
}


// Computes the Laplacian of Psi5 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi5.
//
// Output:  var - the Laplacian of Psi5 for the given particle.
var FermionNumerical::LaplPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    // Using automatic differentiation to compute d^2 psi5 / d x^2  and d^2 psi5 / d y^2.
    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    auto [psi5_xx] = derivativesx(psi5_x, wrt(x(idx)));
    auto [psi5_yy] = derivativesx(psi5_y, wrt(x(idx + 1)));

    return psi5_xx + psi5_yy;
}


// Computes the Laplacian of Psi6 for the specified particle using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of Psi6.
//
// Output:  var - the Laplacian of Psi6 for the given particle.
var FermionNumerical::LaplPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    // Using automatic differentiation to compute d^2 psi6 / d x^2  and d^2 psi6 / d y^2.
    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    auto [psi6_xx] = derivativesx(psi6_x, wrt(x(idx)));
    auto [psi6_yy] = derivativesx(psi6_y, wrt(x(idx + 1)));

    return psi6_xx + psi6_yy;
}


        // Computes the Slater determinant or its gradient component or Laplacian depending on the derivative order.
        //
        // Input:   VectorXvar& x - state vector containing the particle positions;
        //          int particles - the particle set (spin-up or spin-down);
        //          int row_changed - row index of the particle that was changed;
        //          int der_order - order of the derivatives to compute (0 for Psi, 1 for gradient, 2 for Laplacian);
        //          int grad_comp - component of the gradient to compute (if applicable).
        //
        // Output:  double - the computed value based on the requested derivative order.
double FermionNumerical::SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp)
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


// Computes the gradient of the Jastrow factor for a given particle.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of the Jastrow factor.
//
// Output:  VectorXvar - the gradient of the Jastrow factor for the given particle.
VectorXvar FermionNumerical::GradiJastrow(VectorXvar& x, int idx)
{
    var J = Jastrow(x);
    
    // Using automatic differentiation to compute d J / d x  and d J / d y.
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = dJdx;
    grad(1) = dJdy;

    return grad;
}


// Computes the gradient of the Pade-Jastrow factor for a given particle.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the gradient of the Pade-Jastrow factor.
//
// Output:  VectorXvar - the gradient of the Pade-Jastrow factor for the given particle.
VectorXvar FermionNumerical::GradiPadeJastrow(VectorXvar& x, int idx)
{
    var J = PadeJastrow(x);

    // Using automatic differentiation to compute d P / d x  and d P / d y.
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = dJdx;
    grad(1) = dJdy;

    return grad;
}


// Computes the Laplacian of the Jastrow factor for a given particle.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of the Jastrow factor.
//
// Output:  var - the Laplacian of the Jastrow factor for the given particle.
var FermionNumerical::LapliJastrow(VectorXvar& x, int idx)
{
    var J = Jastrow(x);

    // Using automatic differentiation to compute d^2 J / d x^2  and d^2 J / d y^2.
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    auto [d2Jdx2] = derivativesx(dJdx, wrt(x(idx)));
    auto [d2Jdy2] = derivativesx(dJdy, wrt(x(idx + 1)));

    return d2Jdx2 + d2Jdy2;
}


// Computes the Laplacian of the Pade-Jastrow factor for a given particle.
//
// Input:   VectorXvar& x - state vector containing the particle positions;
//          int idx - index of the particle for which to compute the Laplacian of the Pade-Jastrow factor.
//
// Output:  var - the Laplacian of the Pade-Jastrow factor for the given particle.
var FermionNumerical::LapliPadeJastrow(VectorXvar& x, int idx)
{
    var J = PadeJastrow(x);

    // Using automatic differentiation to compute d^2 P / d x^2  and d^2 P / d y^2.
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    auto [d2Jdx2] = derivativesx(dJdx, wrt(x(idx)));
    auto [d2Jdy2] = derivativesx(dJdy, wrt(x(idx + 1)));

    return d2Jdx2 + d2Jdy2;
}


// Computes the Laplacian of PsiT divided by PsiT for all particles in the system using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions.
//
// Output:  double - the computed sum of Laplacians of PsiT divided by PsiT for all particles.
double FermionNumerical::LaplPsiTOverPsiT(VectorXvar& x)
{
    double Psi_up = SD(x, 0, 0, 0, 0);      // Spin-up Slater Determinant.
    double Psi_down = SD(x, 1, 0, 0, 0);    // Spin-down Slater Determinant.

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


// Computes the Laplacian of the Jastrow factor divided by the Jastrow factor for all particles in the system using automatic differentiation.
//
// Input:   VectorXvar& x - state vector containing the particle positions.
//
// Output:  double - the computed sum of Laplacians of Jastrow divided by Jastrow for all particles.
double FermionNumerical::LaplJOverJ(VectorXvar& x)
{
    // If m_mode is 0, compute the Jastrow factor J, otherwise compute the Pade-Jastrow factor P.
    double J = (m_mode == 0) ? val(Jastrow(x)) : val(PadeJastrow(x));

    // This block computes the sum of Laplacians for all particles in the system.
    // If m_mode is 0 (Jastrow ansatz), it computes the Laplacian of the Jastrow factor for each particle.
    // If m_mode is not 0 (Pade-Jastrow ansatz), it computes the Laplacian of the Pade-Jastrow factor for each particle.
    // The loop iterates through all particles, computing the Laplacian for each particle and accumulating the result in the 'sum' variable.
    double sum = 0.0;
    if(m_mode == 0)
    {
        for (int i = 0; i < m_particles; i++)
        {
            int idx = i * 2;
            sum += val(LapliJastrow(x, idx));
        }
    }
    else
    {
        for (int i = 0; i < m_particles; i++)
        {
            int idx = i * 2;
            sum += val(LapliPadeJastrow(x, idx));
        }
    }

    return sum / J;
}


// Computes the gradient of PsiT multiplied by the gradient of the Jastrow factor divided by PsiT for all particles in the system.
//
// Input:   VectorXvar& x - state vector containing the particle positions.
//
// Output:  double - the computed sum of the gradient of PsiT multiplied by the gradient of Jastrow over PsiT for all particles.
double FermionNumerical::GradPsiTGradJOverPsiT(VectorXvar& x)
{

    double Psi_up = SD(x, 0, 0, 0, 0);      // Spin-up Slater Determinant.
    double Psi_down = SD(x, 1, 0, 0, 0);    // Spin-down Slater Determinant.

    // If m_mode is 0, compute the Jastrow factor J, otherwise compute the Pade-Jastrow factor P.
    double J = (m_mode == 0) ? val(Jastrow(x)) : val(PadeJastrow(x));


    // This block computes the sum of gradients for all particles in the system, considering both the spin-up and spin-down components of the wave function.
    // For each particle, the gradient of the Jastrow (or Pade-Jastrow) factor is computed. Based on whether the particle is spin-up or spin-down,
    // the corresponding derivative of the wave function (Psi) is computed and multiplied by the gradient of the Jastrow factor. The results are accumulated in 'sum'.
    // The final result is the sum divided by the product of the spin-up and spin-down wave functions and the Jastrow factor, multiplied by 2.
    double sum = 0.0;
    for (int i = 0; i < m_particles; i++)
    {
        int idx = i * 2;

        // If m_mode is 0, compute the gradient of J (Jastrow factor), otherwise compute the gradient of Pade-Jastrow (PJ) factor.
        VectorXvar gradJ = (m_mode == 0) ? GradiJastrow(x, idx) : GradiPadeJastrow(x, idx);

        double gradJ_x = val(gradJ(0));
        double gradJ_y = val(gradJ(1));

        if (i < m_particles / 2)
        {
            // SD(particles, 0, part_idx, 1, 0) computes the x component of the gradient of the spin-up Slater determinant.
            double dPsi_dx = SD(x, 0, i, 1, 0);

            // SD(particles, 0, part_idx, 1, 1) computes the y component of the gradient of the spin-up Slater determinant.
            double dPsi_dy = SD(x, 0, i, 1, 1);

            sum += dPsi_dx * Psi_down * gradJ_x + dPsi_dy * Psi_down * gradJ_y;
        }
        else
        {
            int j = i - m_particles / 2;

            // SD(particles, 1, part_idx, 1, 0) computes the x component of the gradient of the spin-down Slater determinant.
            double dPsi_dx = SD(x, 1, j, 1, 0);

            // SD(particles, 1, part_idx, 1, 1) computes the y component of the gradient of the spin-down Slater determinant.
            double dPsi_dy = SD(x, 1, j, 1, 1);

            sum += Psi_up * dPsi_dx * gradJ_x + Psi_up * dPsi_dy * gradJ_y;
        }
    }

    double Psi = Psi_up * Psi_down * J;
    return 2.0 * sum / Psi;
}


// Computes the value of the wave function for the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed value of the wave function for the given particle configuration.
double FermionNumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);       // Fill the x vector.

    double Psi_up = SD(x, 0, 0, 0, 0);      // Spin-up Slater Determinant.
    double Psi_down = SD(x, 1, 0, 0, 0);    // Spin-down Slater Determinant.

    // If m_mode is 0, compute the Jastrow factor J, otherwise compute the Pade-Jastrow factor P.
    double J = (m_mode == 0) ? val(Jastrow(x)) : val(PadeJastrow(x));

    return Psi_up * Psi_down * J;
}


// Computes the sum of second derivatives (Laplace operator) of the wave function for the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system.
//
// Output:  double - the computed sum of second derivatives of the wave function for the given particle configuration.
double FermionNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);       // Fill the x vector.

    // Compute all three terms.
    double laplacianPsi1 = LaplPsiTOverPsiT(x);
    double laplacianJ = LaplJOverJ(x);
    double gradgrad = GradPsiTGradJOverPsiT(x);

    return laplacianPsi1 + laplacianJ + gradgrad;
}


// Computes the quantum force for a specific particle in the fermionic system.
//
// Input:   std::vector<std::unique_ptr<class Particle>>& particles - list of particles in the system;
//          int i - index of the particle for which to compute the quantum force.
//
// Output:  std::vector<double> - the quantum force for the specified particle.
std::vector<double> FermionNumerical::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int i)
{
    VectorXvar x = fill_x(particles);   // Fill the x vector.

    int half      = m_particles/2;
    int spinBlock = (i < half ? 0 : 1);     // If 'i' is less than half, it belongs to the spin-up block (spinBlock = 0), otherwise spin-down (spinBlock = 1).
    int localIdx  = (i < half ? i : i - half);  // If 'i' is spin-up, localIdx is the same, otherwise, for spin-down, subtract half of the particles to adjust the index.

    double D_block = SD(x, spinBlock, 0, 0, 0); // Compute the determinant of the Slater determinant for the specified spin block and particle configuration.
    double dlnD_dx = val( SD(x, spinBlock, localIdx, 1, 0) ) / D_block; // Compute the x-component of the gradient of the Slater determinant for the specified particle.
    double dlnD_dy = val( SD(x, spinBlock, localIdx, 1, 1) ) / D_block; // Compute the y-component of the gradient of the Slater determinant for the specified particle.

    VectorXvar gradJ_var;
    if (m_mode == 0)
    {
        gradJ_var = GradiJastrow(x, i*2);       // Compute the gradient of the Jastrow factor for the particle 'i'.
    } else
    {
        gradJ_var = GradiPadeJastrow(x, i*2);   // Compute the gradient of the Pade-Jastrow factor for the particle 'i'.
    }

    // If m_mode is 0, compute the Jastrow factor J, otherwise compute the Pade-Jastrow factor P.
    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);
    double dlnJ_dx = val( gradJ_var(0) ) / val(J);
    double dlnJ_dy = val( gradJ_var(1) ) / val(J);

    return { 2.0 * (dlnD_dx + dlnJ_dx), 2.0 * (dlnD_dy + dlnJ_dy) };
}