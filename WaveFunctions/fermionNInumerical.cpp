#include "fermionNInumerical.h"


FermionNInumerical::FermionNInumerical(double alpha, int n_particles, double omega)
{
    assert(alpha >= 0);

    m_numberOfParameters = 1;

    m_parameters.push_back(alpha);

    m_particles = n_particles;

    m_omega = omega;
    m_sqrt_om = sqrt(omega);
}


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


var FermionNInumerical::psi1i(VectorXvar& x, int idx)
{
    var alpha = m_parameters.back();

    var x_ = x(idx);
    var y_ = x(idx + 1);
    var sum = x_ * x_ + y_ * y_;

    return exp(-alpha * m_omega * sum);
}


var FermionNInumerical::psi2i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx);
}


var FermionNInumerical::psi3i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx + 1);
}

var FermionNInumerical::psi4i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    return psi1 * x(idx) * x(idx + 1);
}


var FermionNInumerical::psi5i(VectorXvar& x, int idx)
{
    var x_ = x(idx);

    return psi1i(x, idx) * (x_ * x_ - 1);
}


var FermionNInumerical::psi6i(VectorXvar& x, int idx)
{
    var y_ = x(idx + 1);

    return psi1i(x, idx) * (y_ * y_ - 1);
}


VectorXvar FermionNInumerical::GradPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi1_x;
    grad(1) = psi1_y;

    return grad;
}


VectorXvar FermionNInumerical::GradPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi2_x;
    grad(1) = psi2_y;

    return grad;
}


VectorXvar FermionNInumerical::GradPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi3_x;
    grad(1) = psi3_y;

    return grad;
}


VectorXvar FermionNInumerical::GradPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi4_x;
    grad(1) = psi4_y;

    return grad;
}


VectorXvar FermionNInumerical::GradPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi5_x;
    grad(1) = psi5_y;

    return grad;
}


VectorXvar FermionNInumerical::GradPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi6_x;
    grad(1) = psi6_y;

    return grad;
}


var FermionNInumerical::LaplPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    auto [psi1_xx] = derivativesx(psi1_x, wrt(x(idx)));
    auto [psi1_yy] = derivativesx(psi1_y, wrt(x(idx + 1)));

    return psi1_xx + psi1_yy;
}


var FermionNInumerical::LaplPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    auto [psi2_xx] = derivativesx(psi2_x, wrt(x(idx)));
    auto [psi2_yy] = derivativesx(psi2_y, wrt(x(idx + 1)));

    return psi2_xx + psi2_yy;
}


var FermionNInumerical::LaplPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    auto [psi3_xx] = derivativesx(psi3_x, wrt(x(idx)));
    auto [psi3_yy] = derivativesx(psi3_y, wrt(x(idx + 1)));

    return psi3_xx + psi3_yy;
}


var FermionNInumerical::LaplPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    auto [psi4_xx] = derivativesx(psi4_x, wrt(x(idx)));
    auto [psi4_yy] = derivativesx(psi4_y, wrt(x(idx + 1)));

    return psi4_xx + psi4_yy;
}


var FermionNInumerical::LaplPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    auto [psi5_xx] = derivativesx(psi5_x, wrt(x(idx)));
    auto [psi5_yy] = derivativesx(psi5_y, wrt(x(idx + 1)));

    return psi5_xx + psi5_yy;
}


var FermionNInumerical::LaplPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    auto [psi6_xx] = derivativesx(psi6_x, wrt(x(idx)));
    auto [psi6_yy] = derivativesx(psi6_y, wrt(x(idx + 1)));

    return psi6_xx + psi6_yy;
}


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


double FermionNInumerical::LaplPsiTOverPsiT(VectorXvar& x)
{
    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);

    double sum = 0.0;
    for (int i = 0; i < m_particles; i++)
    {
        if (i < m_particles / 2)
        {
            sum += SD(x, 0, i, 2, 0) * Psi_down;
        }
        else
        {
            int j = i - m_particles / 2;
            sum += Psi_up * SD(x, 1, j, 2, 0);
        }
    }

    return sum / (Psi_up * Psi_down);
}


double FermionNInumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);

    return Psi_up * Psi_down;
}


double FermionNInumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    return LaplPsiTOverPsiT(x);
}


std::vector<double> FermionNInumerical::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    VectorXvar x = fill_x(particles);
    int size = m_particles / 2;


    double det_up = SD(x, 0, 0, 0, 0);
    double det_down = SD(x, 1, 0, 0, 0);

    std::vector<double> F;

    if(n < m_particles/2)
    {
        double grad_x = SD(x, 0, n, 1, 0);
        double grad_y = SD(x, 0, n, 1, 1);

        F = {2 * grad_x / det_up, 2 * grad_y / det_up};
    }
    else
    {
        double grad_x = SD(x, 1, n - size, 1, 0);
        double grad_y = SD(x, 1, n - size, 1, 1);

        F = {2 * grad_x / det_down, 2 * grad_y / det_down};
    }

    return F;
}