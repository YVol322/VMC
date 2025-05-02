#include "fermionnumerical.h"


FermionNumerical::FermionNumerical(double alpha, std::vector<double>beta, int mode, int n_particles)
{
    assert(alpha >= 0);

    int n_betas = beta.size();
    m_numberOfParameters = n_betas + 1;

    m_parameters.reserve(m_numberOfParameters);
    m_parameters.insert(m_parameters.end(), beta.begin(), beta.end());
    m_parameters.push_back(alpha);

    m_particles = n_particles;
    m_mode = mode;
}


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


int FermionNumerical::BetaIndex(int i, int j)
{
    int n_particles = m_particles;

    int ii = std::min(i, j);
    int jj = std::max(i, j);

    return ii * (2 * n_particles - ii - 1) / 2 + (jj - ii - 1);
}


var FermionNumerical::a_ij(int i, int j)
{
    int n_particles = m_particles;

    double aij;
    
    if ((i < n_particles / 2 && j >= n_particles / 2) || (i >= n_particles / 2 && j < n_particles / 2))
    {
        aij = 1;
    }
    else
    {
        aij = 1.0 / 3.0;
    }

    return aij;
}


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
            var rij = sqrt(sum2);

            int idx = BetaIndex(i,j);
            var beta_ij = m_parameters[idx];

            sum1 += beta_ij * rij;
        }
    }

    return exp(sum1);
}


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
            var rij = sqrt(sum2);

            var aij = a_ij(i, j);

            sum1 += aij * rij / (1 + beta * rij);
        }
    }

    return exp(sum1);
}


var FermionNumerical::psi1i(VectorXvar& x, int idx)
{
    var alpha = m_parameters.back();

    var x_ = x(idx);
    var y_ = x(idx + 1);
    var sum = x_ * x_ + y_ * y_;

    return exp(-alpha * sum);
}


var FermionNumerical::psi2i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx);
}


var FermionNumerical::psi3i(VectorXvar& x, int idx)
{
    return psi1i(x, idx) * x(idx + 1);
}

var FermionNumerical::psi4i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    return psi1 * x(idx) * x(idx + 1);
}


var FermionNumerical::psi5i(VectorXvar& x, int idx)
{
    var x_ = x(idx);

    return psi1i(x, idx) * (x_ * x_ - 1);
}


var FermionNumerical::psi6i(VectorXvar& x, int idx)
{
    var y_ = x(idx + 1);

    return psi1i(x, idx) * (y_ * y_ - 1);
}


VectorXvar FermionNumerical::GradPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi1_x;
    grad(1) = psi1_y;

    return grad;
}

VectorXvar FermionNumerical::GradPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi2_x;
    grad(1) = psi2_y;

    return grad;
}


VectorXvar FermionNumerical::GradPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi3_x;
    grad(1) = psi3_y;

    return grad;
}


VectorXvar FermionNumerical::GradPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi4_x;
    grad(1) = psi4_y;

    return grad;
}


VectorXvar FermionNumerical::GradPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi5_x;
    grad(1) = psi5_y;

    return grad;
}


VectorXvar FermionNumerical::GradPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = psi6_x;
    grad(1) = psi6_y;

    return grad;
}


var FermionNumerical::LaplPsi1i(VectorXvar& x, int idx)
{
    var psi1 = psi1i(x, idx);

    auto [psi1_x, psi1_y] = derivativesx(psi1, wrt(x(idx), x(idx + 1)));

    auto [psi1_xx] = derivativesx(psi1_x, wrt(x(idx)));
    auto [psi1_yy] = derivativesx(psi1_y, wrt(x(idx + 1)));

    return psi1_xx + psi1_yy;
}


var FermionNumerical::LaplPsi2i(VectorXvar& x, int idx)
{
    var psi2 = psi2i(x, idx);

    auto [psi2_x, psi2_y] = derivativesx(psi2, wrt(x(idx), x(idx + 1)));

    auto [psi2_xx] = derivativesx(psi2_x, wrt(x(idx)));
    auto [psi2_yy] = derivativesx(psi2_y, wrt(x(idx + 1)));

    return psi2_xx + psi2_yy;
}


var FermionNumerical::LaplPsi3i(VectorXvar& x, int idx)
{
    var psi3 = psi3i(x, idx);

    auto [psi3_x, psi3_y] = derivativesx(psi3, wrt(x(idx), x(idx + 1)));

    auto [psi3_xx] = derivativesx(psi3_x, wrt(x(idx)));
    auto [psi3_yy] = derivativesx(psi3_y, wrt(x(idx + 1)));

    return psi3_xx + psi3_yy;
}


var FermionNumerical::LaplPsi4i(VectorXvar& x, int idx)
{
    var psi4 = psi4i(x, idx);

    auto [psi4_x, psi4_y] = derivativesx(psi4, wrt(x(idx), x(idx + 1)));

    auto [psi4_xx] = derivativesx(psi4_x, wrt(x(idx)));
    auto [psi4_yy] = derivativesx(psi4_y, wrt(x(idx + 1)));

    return psi4_xx + psi4_yy;
}


var FermionNumerical::LaplPsi5i(VectorXvar& x, int idx)
{
    var psi5 = psi5i(x, idx);

    auto [psi5_x, psi5_y] = derivativesx(psi5, wrt(x(idx), x(idx + 1)));

    auto [psi5_xx] = derivativesx(psi5_x, wrt(x(idx)));
    auto [psi5_yy] = derivativesx(psi5_y, wrt(x(idx + 1)));

    return psi5_xx + psi5_yy;
}


var FermionNumerical::LaplPsi6i(VectorXvar& x, int idx)
{
    var psi6 = psi6i(x, idx);

    auto [psi6_x, psi6_y] = derivativesx(psi6, wrt(x(idx), x(idx + 1)));

    auto [psi6_xx] = derivativesx(psi6_x, wrt(x(idx)));
    auto [psi6_yy] = derivativesx(psi6_y, wrt(x(idx + 1)));

    return psi6_xx + psi6_yy;
}


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


VectorXvar FermionNumerical::GradiJastrow(VectorXvar& x, int idx)
{
    var J = Jastrow(x);
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = dJdx;
    grad(1) = dJdy;

    return grad;
}


VectorXvar FermionNumerical::GradiPadeJastrow(VectorXvar& x, int idx)
{
    var J = PadeJastrow(x);
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    VectorXvar grad(2);
    grad(0) = dJdx;
    grad(1) = dJdy;

    return grad;
}


var FermionNumerical::LapliJastrow(VectorXvar& x, int idx)
{
    var J = Jastrow(x);

    // First derivatives
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    auto [d2Jdx2] = derivativesx(dJdx, wrt(x(idx)));
    auto [d2Jdy2] = derivativesx(dJdy, wrt(x(idx + 1)));

    return d2Jdx2 + d2Jdy2;
}


var FermionNumerical::LapliPadeJastrow(VectorXvar& x, int idx)
{
    var J = PadeJastrow(x);

    // First derivatives
    auto [dJdx, dJdy] = derivativesx(J, wrt(x(idx), x(idx + 1)));

    auto [d2Jdx2] = derivativesx(dJdx, wrt(x(idx)));
    auto [d2Jdy2] = derivativesx(dJdy, wrt(x(idx + 1)));

    return d2Jdx2 + d2Jdy2;
}


double FermionNumerical::LaplacianPsi1OverPsi1(VectorXvar& x)
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


double FermionNumerical::LaplacianJOverJ(VectorXvar& x)
{
    double J = (m_mode == 0) ? val(Jastrow(x)) : val(PadeJastrow(x));

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


double FermionNumerical::GradPsi1GradJOverPsi(VectorXvar& x)
{
    double sum = 0.0;

    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);
    double J = (m_mode == 0) ? val(Jastrow(x)) : val(PadeJastrow(x));

    for (int i = 0; i < m_particles; i++)
    {
        int idx = i * 2;

        VectorXvar gradJ = (m_mode == 0) ? GradiJastrow(x, idx) : GradiPadeJastrow(x, idx);

        double gradJ_x = val(gradJ(0));
        double gradJ_y = val(gradJ(1));

        if (i < m_particles / 2)
        {
            double dPsi_dx = SD(x, 0, i, 1, 0);
            double dPsi_dy = SD(x, 0, i, 1, 1);

            sum += dPsi_dx * Psi_down * gradJ_x + dPsi_dy * Psi_down * gradJ_y;
        }
        else
        {
            int j = i - m_particles / 2;
            double dPsi_dx = SD(x, 1, j, 1, 0);
            double dPsi_dy = SD(x, 1, j, 1, 1);

            sum += Psi_up * dPsi_dx * gradJ_x + Psi_up * dPsi_dy * gradJ_y;
        }
    }

    double Psi = Psi_up * Psi_down * J;
    return 2.0 * sum / Psi;
}


double FermionNumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    double Psi_up = SD(x, 0, 0, 0, 0);
    double Psi_down = SD(x, 1, 0, 0, 0);
    double J = (m_mode == 0) ? val(Jastrow(x)) : val(PadeJastrow(x));

    return Psi_up * Psi_down * J;
}


double FermionNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    double laplacianPsi1 = LaplacianPsi1OverPsi1(x);
    double laplacianJ = LaplacianJOverJ(x);
    double gradgrad = GradPsi1GradJOverPsi(x);

    return laplacianPsi1 + laplacianJ + gradgrad;
}


double FermionNumerical::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    Particle& particle_i = *(particles[part_idx]);
    int n_dimensions = particle_i.getNumberOfDimensions();
    double coordinate;
    double r2 = 0;

    for(int i = 0; i < n_dimensions; i++)
    {
        coordinate = particle_i.getPosition()[i];
        r2 += coordinate * coordinate;
    }

    return r2;
}

double FermionNumerical::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
{
    Particle& pi = *(particles[i]);
    Particle& pj = *(particles[j]);

    int n_dimensions = pi.getNumberOfDimensions();

    double rij = 0;
    for(int d = 0; d < n_dimensions; d++)
    {
        rij += (pi.getPosition()[d] - pj.getPosition()[d]) * (pi.getPosition()[d] - pj.getPosition()[d]);
    }
    
    return sqrt(rij);
}