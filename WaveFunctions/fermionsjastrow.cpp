#include <cmath>
//#include <armadillo>
#include <Eigen/Dense>
#include <iostream>


#include "fermionsjastrow.h"

FermionsJastrow::FermionsJastrow(double alpha, std::vector<double> beta, int mode, int n_particles)
{
    assert(alpha >= 0);
    int n_betas = beta.size();
    m_numberOfParameters = n_betas + 1;
    m_parameters.reserve(m_numberOfParameters);

    m_particles = n_particles;

    for(int i = 0; i < n_betas; i++)
    {
        m_parameters.push_back(beta.at(i));
    }

    m_parameters.push_back(alpha);

    m_mode = mode;
}


double FermionsJastrow::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int N = m_particles;

    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            double rij = r_ij(particles, i, j);

            int idx = i * (2 * N - i - 1) / 2 + (j - i - 1);
            double beta_ij = m_parameters.at(idx);

            sum += beta_ij * rij;
        }
    }

    return std::exp(sum);
}


double FermionsJastrow::PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int N = m_particles;
    double aij = 0;

    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            if ((i < N/2 && j >= N/2) || (i >= N/2 && j < N/2))
            {
                aij = 1;
            }
            else
            {
                aij = 1.0/3.0;
            }
            
            double rij = r_ij(particles, i, j);
            double beta = m_parameters.at(0);

            sum += aij * rij / (1 + beta * rij);
        }
    }

    return std::exp(sum);
}


std::vector<double> FermionsJastrow::GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    int n_particles = particles.size();
    int n_dimensions = particles[0]->getNumberOfDimensions();

    Particle particle_i = *(particles.at(part_inx));

    std::vector<double> grad(n_dimensions, 0.0);

    for(int l = 0; l < n_particles; l++)
    {
        if(l == part_inx) continue;

        Particle particle_l = *(particles.at(l));
        double ril = r_ij(particles, part_inx, l);

        int i = std::min(part_inx, static_cast<double>(l));
        int j = std::max(part_inx, static_cast<double>(l));
        int idx = i * (2 * n_particles - i - 1) / 2 + (j - i - 1);
        double beta_il = m_parameters.at(idx);

        for(int k = 0; k < n_dimensions; k++)
        {
            double coord_i = particle_i.getPosition().at(k);
            double coord_l = particle_l.getPosition().at(k);

            grad.at(k) += beta_il * (coord_i - coord_l) / ril;
        }
    }

    return grad;
}


std::vector<double> FermionsJastrow::GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    int n_particles = particles.size();
    int n_dimensions = particles[0]->getNumberOfDimensions();

    Particle particle_i = *(particles.at(part_inx));

    std::vector<double> grad(n_dimensions, 0.0);

    double aij = 0;

    for(int l = 0; l < n_particles; l++)
    {
        if(l == part_inx) continue;

        Particle particle_l = *(particles.at(l));
        double ril = r_ij(particles, part_inx, l);

        double beta = m_parameters.at(0);

        if ((part_inx < n_particles/2 && l >= n_particles/2) || (part_inx >= n_particles/2 && l < n_particles/2))
        {
            aij = 1;
        }
        else
        {
            aij = 1.0/3.0;
        }

        for(int k = 0; k < n_dimensions; k++)
        {
            double coord_i = particle_i.getPosition().at(k);
            double coord_l = particle_l.getPosition().at(k);

            grad.at(k) += aij * (coord_i - coord_l) / (ril * (1 + beta * ril) * (1 + beta * ril));
        }
    }

    return grad;
}


double FermionsJastrow::Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));
    int n_dimensions = particle_i.getNumberOfDimensions();
    double r2 = 0;
    double cord, psi1;

    for(int i = 0; i < n_dimensions; i++)
    {
        cord = (particle_i.getPosition()).at(i);
        r2 += cord * cord;
    }

    r2 *= -alpha;
    psi1 = exp(r2);

    return psi1;
}


double FermionsJastrow::Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    Particle particle_i = *(particles.at(part_inx));
    double x = particle_i.getPosition().at(0);
    double psi2 = x * Psi1(particles, part_inx);

    return psi2;
}

double FermionsJastrow::Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    Particle particle_i = *(particles.at(part_inx));
    double y = particle_i.getPosition().at(1);
    double psi3 = y * Psi1(particles, part_inx);

    return psi3;
}


double FermionsJastrow::Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    Particle particle_i = *(particles.at(part_inx));
    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi4 = x * y * Psi1(particles, part_inx);

    return psi4;
}


double FermionsJastrow::Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    Particle particle_i = *(particles.at(part_inx));
    double x = particle_i.getPosition().at(0);
    double psi5 = (x * x - 1) * Psi1(particles, part_inx);

    return psi5;
}


double FermionsJastrow::Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    Particle particle_i = *(particles.at(part_inx));
    double y = particle_i.getPosition().at(1);
    double psi6 = (y * y - 1) * Psi1(particles, part_inx);

    return psi6;
}

std::vector<double> FermionsJastrow::GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double dpsi1dx = -2 * alpha * x * psi1;
    double dpsi1dy = -2 * alpha * y * psi1;

    std::vector<double> grad_psi1(2);

    grad_psi1.at(0) = dpsi1dx;
    grad_psi1.at(1) = dpsi1dy;

    return grad_psi1;
}

std::vector<double> FermionsJastrow::GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double dpsi2dx = (1 - 2 * alpha * x * x) * psi1;
    double dpsi2dy = -2 * alpha * y * x * psi1;

    std::vector<double> grad_psi2(2);

    grad_psi2.at(0) = dpsi2dx;
    grad_psi2.at(1) = dpsi2dy;

    return grad_psi2;
}

std::vector<double> FermionsJastrow::GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double dpsi3dx = -2 * alpha * x * y * psi1;
    double dpsi3dy = (1 - 2 * alpha * y * y) * psi1;

    std::vector<double> grad_psi3(2);

    grad_psi3.at(0) = dpsi3dx;
    grad_psi3.at(1) = dpsi3dy;

    return grad_psi3;
}


std::vector<double> FermionsJastrow::GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double dpsi4dx = y * psi1 * (1 - 2 * alpha * x * x);
    double dpsi4dy = x * psi1 * (1 - 2 * alpha * y * y);

    std::vector<double> grad_psi4(2);

    grad_psi4.at(0) = dpsi4dx;
    grad_psi4.at(1) = dpsi4dy;

    return grad_psi4;
}


std::vector<double> FermionsJastrow::GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double dpsi5dx = 2 * x * psi1 * (1 - alpha * x * x + alpha);
    double dpsi5dy = -2 * alpha * y * (x * x - 1) * psi1;

    std::vector<double> grad_psi5(2);

    grad_psi5.at(0) = dpsi5dx;
    grad_psi5.at(1) = dpsi5dy;

    return grad_psi5;
}


std::vector<double> FermionsJastrow::GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double dpsi6dx = -2 * alpha * x * (y * y - 1) * psi1;
    double dpsi6dy = 2 * y * psi1 * (1 - alpha * y * y + alpha);

    std::vector<double> grad_psi6(2);

    grad_psi6.at(0) = dpsi6dx;
    grad_psi6.at(1) = dpsi6dy;

    return grad_psi6;
}


double FermionsJastrow::LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi1 = (-4 * alpha + 4 * alpha * alpha * (x * x + y * y) ) * psi1;

    return laplacian_psi1;
}

double FermionsJastrow::LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi2 = (-8 * alpha * x + 4 * alpha * alpha * x * (x * x + y * y) ) * psi1;

    return laplacian_psi2;
}

double FermionsJastrow::LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi3 = (-8 * alpha * y + 4 * alpha * alpha * y * (x * x + y * y) ) * psi1;

    return laplacian_psi3;
}


double FermionsJastrow::LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double r2 = x * x + y * y;
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi4 = x * y * psi1 * (4 * alpha * alpha * r2 - 12 * alpha);

    return laplacian_psi4;
}


double FermionsJastrow::LapliPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double x2 = x * x;
    double y2 = y * y;
    double common = x2 - 1;

    double laplacian_psi5 = 2 * psi1 * (
        -4 * alpha * x2
        + alpha * common * (2 * alpha * x2 - 1)
        + alpha * common * (2 * alpha * y2 - 1)
        + 1
    );

    return laplacian_psi5;
}


double FermionsJastrow::LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double x2 = x * x;
    double y2 = y * y;
    double common = y2 - 1;

    double laplacian_psi6 = 2 * psi1 * (
        -4 * alpha * y2
        + alpha * common * (2 * alpha * x2 - 1)
        + alpha * common * (2 * alpha * y2 - 1)
        + 1
    );

    return laplacian_psi6;
}

double FermionsJastrow::SD(std::vector<std::unique_ptr<class Particle>>& particles, int particles_set, int row_changed, int der_order, int grad_comp)
{
    int N = m_particles;
    double det = 0;

    if(N/2 == 3)
    {
        Eigen::MatrixXd A(N/2, N/2);
        for(int i = 0; i < N/2; i++)
        {
            if (i == row_changed)
            {
                if (der_order == 0)
                {
                    A(i, 0) = Psi1(particles, i + N/2 * particles_set);
                    A(i, 1) = Psi2(particles, i + N/2 * particles_set);
                    A(i, 2) = Psi3(particles, i + N/2 * particles_set);
                }
                else if (der_order == 1)
                {
                    A(i, 0) = GradiPsi1(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 1) = GradiPsi2(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 2) = GradiPsi3(particles, i + N/2 * particles_set).at(grad_comp);
                }
                else if (der_order == 2)
                {
                    A(i, 0) = LapliPsi1(particles, i + N/2 * particles_set);
                    A(i, 1) = LapliPsi2(particles, i + N/2 * particles_set);
                    A(i, 2) = LapliPsi3(particles, i + N/2 * particles_set);
                }
            }
            else
            {
                A(i, 0) = Psi1(particles, i + N/2 * particles_set);
                A(i, 1) = Psi2(particles, i + N/2 * particles_set);
                A(i, 2) = Psi3(particles, i + N/2 * particles_set);
            }
        }
        det = A.determinant();
    }
    else if(N/2 == 1)
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
    if(N/2 == 6)
    {
        Eigen::MatrixXd A(N/2, N/2);
        for(int i = 0; i < N/2; i++)
        {
            if (i == row_changed)
            {
                if (der_order == 0)
                {
                    A(i, 0) = Psi1(particles, i + N/2 * particles_set);
                    A(i, 1) = Psi2(particles, i + N/2 * particles_set);
                    A(i, 2) = Psi3(particles, i + N/2 * particles_set);
                    A(i, 3) = Psi4(particles, i + N/2 * particles_set);
                    A(i, 4) = Psi5(particles, i + N/2 * particles_set);
                    A(i, 5) = Psi6(particles, i + N/2 * particles_set);
                }
                else if (der_order == 1)
                {
                    A(i, 0) = GradiPsi1(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 1) = GradiPsi2(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 2) = GradiPsi3(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 3) = GradiPsi4(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 4) = GradiPsi5(particles, i + N/2 * particles_set).at(grad_comp);
                    A(i, 5) = GradiPsi6(particles, i + N/2 * particles_set).at(grad_comp);

                }
                else if (der_order == 2)
                {
                    A(i, 0) = LapliPsi1(particles, i + N/2 * particles_set);
                    A(i, 1) = LapliPsi2(particles, i + N/2 * particles_set);
                    A(i, 2) = LapliPsi3(particles, i + N/2 * particles_set);
                    A(i, 3) = LapliPsi4(particles, i + N/2 * particles_set);
                    A(i, 4) = LapliPsi5(particles, i + N/2 * particles_set);
                    A(i, 5) = LapliPsi6(particles, i + N/2 * particles_set);
                }
            }
            else
            {
                A(i, 0) = Psi1(particles, i + N/2 * particles_set);
                A(i, 1) = Psi2(particles, i + N/2 * particles_set);
                A(i, 2) = Psi3(particles, i + N/2 * particles_set);
                A(i, 3) = Psi4(particles, i + N/2 * particles_set);
                A(i, 4) = Psi5(particles, i + N/2 * particles_set);
                A(i, 5) = Psi6(particles, i + N/2 * particles_set);
            }
        }
        det = A.determinant();
    }

    return det;
}


double FermionsJastrow::LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0]->getNumberOfDimensions();

    std::vector<double> grad_i = GradiJOverJ(particles, part_inx);
    double sum1 = 0.0;
    for (int m = 0; m < n_dimensions; ++m)
    {
        sum1 += grad_i.at(m) * grad_i.at(m);
    }

    double sum2 = 0.0;
    for (int l = 0; l < n_particles; ++l)
    {
        if (l == part_inx) continue;

        int i = std::min(part_inx, static_cast<double>(l));
        int j = std::max(part_inx, static_cast<double>(l));
        int idx = i * (2 * n_particles - i - 1) / 2 + (j - i - 1);
        double beta_il = m_parameters.at(idx);

        double ril = r_ij(particles, part_inx, l);

        sum2 += beta_il / ril;
    }

    return sum1 + sum2;
}


double FermionsJastrow::LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0]->getNumberOfDimensions();

    std::vector<double> grad_i = GradiPJOverPJ(particles, part_inx);
    double sum1 = 0.0;
    for (int m = 0; m < n_dimensions; ++m)
    {
        sum1 += grad_i.at(m) * grad_i.at(m);
    }

    double sum2 = 0.0;
    double aij = 0.0;
    double beta = m_parameters.at(0);

    for (int l = 0; l < n_particles; ++l)
    {
        if (l == part_inx) continue;

        if ((part_inx < n_particles/2 && l >= n_particles/2) || (part_inx >= n_particles/2 && l < n_particles/2))
        {
            aij = 1;
        }
        else
        {
            aij = 1.0/3.0;
        }

        double ril = r_ij(particles, part_inx, l);
        double t = 1 + beta * ril;

        sum2 += aij * ( 1/(ril * t * t) - 2*beta/(t*t*t) );
    }

    return sum1 + sum2;
}



double FermionsJastrow::LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles)
{
    int n_particles = m_particles;

    double det_up = SD(particles, 0, 0, 0, 0);
    double det_down = SD(particles, 1, 0, 0, 0);

    double sum = 0.0;

    for(int i = 0; i < n_particles/2; i++)
    {
        sum += SD(particles, 0, i, 2, 0) / det_up;
        sum += SD(particles, 1, i, 2, 0) / det_down;
    }

    return sum;
}


double FermionsJastrow::GradiPsi1GradiJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0]->getNumberOfDimensions();

    std::vector<double> grad_J;
    if(m_mode == 0) grad_J = GradiJOverJ(particles, part_inx);
    else grad_J = GradiPJOverPJ(particles, part_inx);


    std::vector<double> grad_psi(n_dimensions);

    double det_up = SD(particles, 0, 0, 0, 0);
    double det_down = SD(particles, 1, 0, 0, 0);

    double sum = 0.0;

    if (part_inx < n_particles/2)
    {
        for (int i = 0; i < n_dimensions; ++i)
        {
            grad_psi.at(i) = SD(particles, 0, part_inx, 1, i) / det_up;
            sum += grad_psi.at(i) * grad_J.at(i);
        }
    }
    else
    {
        int local_index = part_inx - n_particles/2;
        for (int i = 0; i < n_dimensions; ++i)
        {
            grad_psi.at(i) = SD(particles, 1, local_index, 1, i) / det_down;
            sum += grad_psi.at(i) * grad_J.at(i);
        }
    }

    return 2.0 * sum;
}



double FermionsJastrow::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double det_up = SD(particles, 0, 0, 0, 0);
    double det_down = SD(particles, 1, 0, 0, 0);
    double J;

    if(m_mode == 0) J = Jastrow(particles);
    else J = PadeJastrow(particles);

    return det_up * det_down * J;
}

double FermionsJastrow::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    int n_particles = particles.size();

    double LaplPsi = LaplPsi1OverPsi1(particles);
    double LaplJ = 0;
    double GradGrad = 0;

    if(m_mode == 0)
    {
        for(int i = 0; i < n_particles; i++)
        {
            LaplJ += LapliJOverJ(particles, i);
            GradGrad += GradiPsi1GradiJOverPsi(particles, i);
        }
    }
    else
    {
        for(int i = 0; i < n_particles; i++)
        {
            LaplJ += LapliPJOverPJ(particles, i);
            GradGrad += GradiPsi1GradiJOverPsi(particles, i);
        }
    }

    return LaplPsi + LaplJ + GradGrad;
}


double FermionsJastrow::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
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

double FermionsJastrow::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
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