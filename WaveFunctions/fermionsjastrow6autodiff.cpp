#include <cmath>
//#include <armadillo>
#include <Eigen/Dense>
#include <iostream>


#include "fermionsjastrow6autodiff.h"

FermionsJastrow6Autodiff::FermionsJastrow6Autodiff(double alpha, std::vector<double>beta)
{
    assert(alpha >= 0);
    m_numberOfParameters = 16;
    m_parameters.reserve(m_numberOfParameters);
    for(int i = 0; i < 15; i++)
    {
        m_parameters.push_back(beta.at(i));
    }

    m_parameters.push_back(alpha);
}

VectorXvar FermionsJastrow6Autodiff::fill_x(std::vector<std::unique_ptr<class Particle>>& particles)
{
    int n_particles = particles.size();
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    VectorXvar x(n_particles * n_dimensions);
    for(int i = 0; i < n_particles; i++)
    {
        for(int j = 0; j < n_dimensions; j++)
        {
            x(i * n_dimensions + j) = particles[i]->getPosition().at(j);
        }
    }
    return x;
}

var FermionsJastrow6Autodiff::psi1i(VectorXvar& x, int idx)
{
    var alpha = m_parameters.back();

    var sum = x(idx) * x(idx) + x(idx + 1) * x(idx + 1);

    sum *= -alpha;

    return exp(sum);
}


var FermionsJastrow6Autodiff::psi2i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    return psi1 * x(idx);
}


var FermionsJastrow6Autodiff::psi3i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    return psi1 * x(idx + 1);
}


VectorXvar FermionsJastrow6Autodiff::GradPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi1_x;
    grad(1) = psi1_y;

    return grad;
}

VectorXvar FermionsJastrow6Autodiff::GradPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi2_x;
    grad(1) = psi2_y;

    return grad;
}


VectorXvar FermionsJastrow6Autodiff::GradPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi3_x;
    grad(1) = psi3_y;

    return grad;
}


var FermionsJastrow6Autodiff::LaplPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    auto [psi1_xx] = derivativesx(psi1_x, wrt(x(idx)));
    auto [psi1_yy] = derivativesx(psi1_y, wrt(x(idx + 1)));

    return psi1_xx + psi1_yy;
}


var FermionsJastrow6Autodiff::LaplPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    auto [psi2_xx] = derivativesx(psi2_x, wrt(x(idx)));
    auto [psi2_yy] = derivativesx(psi2_y, wrt(x(idx + 1)));

    return psi2_xx + psi2_yy;
}


var FermionsJastrow6Autodiff::LaplPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    auto [psi3_xx] = derivativesx(psi3_x, wrt(x(idx)));
    auto [psi3_yy] = derivativesx(psi3_y, wrt(x(idx + 1)));

    return psi3_xx + psi3_yy;
}


var FermionsJastrow6Autodiff::Jastrow(const VectorXvar& x)
{
    int N = x.size() / 2;
    var sum = 0.0;

    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            var dx = x(2 * i) - x(2 * j);
            var dy = x(2 * i + 1) - x(2 * j + 1);
            var rij = sqrt(dx * dx + dy * dy);

            int idx = i * (2 * N - i - 1) / 2 + (j - i - 1);
            var beta_ij = m_parameters.at(idx);

            sum += beta_ij * rij;
        }
    }

    return exp(sum);
}

VectorXvar FermionsJastrow6Autodiff::GradiJastrow(VectorXvar& x, int idx)
{
    auto wrapped_Jastrow = [&](const VectorXvar& x_) {
        return Jastrow(x_);
    };

    var J = wrapped_Jastrow(x);
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = dJdx;
    grad(1) = dJdy;

    return grad;
}

var FermionsJastrow6Autodiff::LapliJastrow(VectorXvar& x, int idx)
{
    // Wrap Jastrow in a lambda for autodiff
    auto wrapped_Jastrow = [&](const VectorXvar& x_) {
        return Jastrow(x_);
    };

    var J = wrapped_Jastrow(x);

    // First derivatives
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    // Second derivatives (Laplacian = d²J/dx² + d²J/dy²)
    auto [d2Jdx2] = derivativesx(dJdx, wrt(x(idx)));
    auto [d2Jdy2] = derivativesx(dJdy, wrt(x(idx + 1)));

    return d2Jdx2 + d2Jdy2;
}

double FermionsJastrow6Autodiff::SD(VectorXvar& x, int particles, int row_changed, int der_order, int grad_comp)
{
    MatrixXvar A(3, 3);

    for(int i = 0; i < 3; i++)
    {
        int idx = i * 2 + 6 * particles;

        if (i == row_changed)
        {
            if (der_order == 0)
            {
                A(i, 0) = psi1i(x, idx);
                A(i, 1) = psi2i(x, idx);
                A(i, 2) = psi3i(x, idx);
            }
            else if (der_order == 1)
            {
                VectorXvar grad1 = GradPsi1i(x, idx);
                VectorXvar grad2 = GradPsi2i(x, idx);
                VectorXvar grad3 = GradPsi3i(x, idx);

                A(i, 0) = grad1(grad_comp);
                A(i, 1) = grad2(grad_comp);
                A(i, 2) = grad3(grad_comp);
            }
            else if (der_order == 2)
            {
                A(i, 0) = LaplPsi1i(x, idx);
                A(i, 1) = LaplPsi2i(x, idx);
                A(i, 2) = LaplPsi3i(x, idx);
            }
        }
        else
        {
            A(i, 0) = psi1i(x, idx);
            A(i, 1) = psi2i(x, idx);
            A(i, 2) = psi3i(x, idx);
        }
    }

    var det = 
      A(0,0) * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    - A(0,1) * (A(1,0)*A(2,2) - A(1,2)*A(2,0))
    + A(0,2) * (A(1,0)*A(2,1) - A(1,1)*A(2,0));

    return val(det);
}


double FermionsJastrow6Autodiff::GradPsi1GradJOverPsi(VectorXvar& x)
{
    double sum = 0.0;

    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);
    double J = val(Jastrow(x));

    for (int i = 0; i < 6; ++i)
    {
        int idx = i * 2;

        VectorXvar gradJ = GradiJastrow(x, idx);
        double gradJ_x = val(gradJ(0));
        double gradJ_y = val(gradJ(1));

        if (i < 3)
        {
            double dPsi_dx = SD(x, 0, i, 1, 0);
            double dPsi_dy = SD(x, 0, i, 1, 1);

            sum += dPsi_dx * Psi_down * gradJ_x + dPsi_dy * Psi_down * gradJ_y;
        }
        else
        {
            int j = i - 3;
            double dPsi_dx = SD(x, 1, j, 1, 0);
            double dPsi_dy = SD(x, 1, j, 1, 1);

            sum += Psi_up * dPsi_dx * gradJ_x + Psi_up * dPsi_dy * gradJ_y;
        }
    }

    double Psi = Psi_up * Psi_down * J;
    return 2.0 * sum / Psi;
}


double FermionsJastrow6Autodiff::LaplacianPsi1OverPsi1(VectorXvar& x)
{
    double sum = 0.0;

    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);

    for (int i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            sum += SD(x, 0, i, 2, 0) * Psi_down;
        }
        else
        {
            int j = i - 3;
            sum += Psi_up * SD(x, 1, j, 2, 0);
        }
    }

    return sum / (Psi_up * Psi_down);
}

double FermionsJastrow6Autodiff::LaplacianJOverJ(VectorXvar& x)
{
    double J = val(Jastrow(x));
    double sum = 0.0;

    for (int i = 0; i < 6; ++i)
    {
        int idx = i * 2;
        sum += val(LapliJastrow(x, idx));
    }

    return sum / J;
}


double FermionsJastrow6Autodiff::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);
    double J = val(Jastrow(x));

    return Psi_up * Psi_down * J;
}



double FermionsJastrow6Autodiff::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    double laplacianPsi1 = LaplacianPsi1OverPsi1(x);
    double laplacianJ = LaplacianJOverJ(x);
    double crossTerm = GradPsi1GradJOverPsi(x);

    double result = laplacianPsi1 + laplacianJ + crossTerm;

    return result;
}

double FermionsJastrow6Autodiff::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
{
    Particle particle_i = *(particles.at(part_index));
    int n_dimensions = particle_i.getNumberOfDimensions();
    double coordinate;
    double r2 = 0;

    for(int i = 0; i < n_dimensions; i++)
    {
        coordinate = particle_i.getPosition().at(i);
        r2 += coordinate * coordinate;
    }

    return r2;
}

double FermionsJastrow6Autodiff::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
{
    Particle& pi = *(particles.at(i));
    Particle& pj = *(particles.at(j));

    double x_i = pi.getPosition().at(0);
    double y_i = pi.getPosition().at(1);
    double x_j = pj.getPosition().at(0);
    double y_j = pj.getPosition().at(1);

    double r_ij = sqrt((x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j));

    return r_ij;
}