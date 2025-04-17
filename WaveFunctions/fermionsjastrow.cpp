#include <cmath>
//#include <armadillo>
#include <Eigen/Dense>
#include <iostream>


#include "fermionsjastrow.h"

FermionsJastrow::FermionsJastrow(double alpha, std::vector<double> beta)
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

double FermionsJastrow::Psi1DoubleDer(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi1 = (-4 * alpha + 4 * alpha * alpha * (x * x + y * y) ) * psi1;

    return laplacian_psi1;
}

double FermionsJastrow::Psi2DoubleDer(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi2 = (-8 * alpha * x + 4 * alpha * alpha * x * (x * x + y * y) ) * psi1;

    return laplacian_psi2;
}

double FermionsJastrow::Psi3DoubleDer(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
{
    double alpha = m_parameters.back();
    Particle particle_i = *(particles.at(part_inx));

    double x = particle_i.getPosition().at(0);
    double y = particle_i.getPosition().at(1);
    double psi1 = Psi1(particles, part_inx);

    double laplacian_psi3 = (-8 * alpha * y + 4 * alpha * alpha * y * (x * x + y * y) ) * psi1;

    return laplacian_psi3;
}

double FermionsJastrow::SD(std::vector<std::unique_ptr<class Particle>>& particles, double row_indx, double parts, int der_order, int coord)
{
    Eigen::MatrixXd A(3, 3);

    for(int i = 0; i < 3; i++)
    {
        if(i == row_indx)
        {
            if(der_order == 2)
            {
                A(i, 0) = Psi1DoubleDer(particles, i + 3 * parts);
                A(i, 1) = Psi2DoubleDer(particles, i + 3 * parts);
                A(i, 2) = Psi3DoubleDer(particles, i + 3 * parts);
            }
            else if(der_order == 1)
            {
                if(coord == 1)
                {
                    A(i, 0) = Psi1Der(particles, i + 3 * parts).at(0);
                    A(i, 1) = Psi2Der(particles, i + 3 * parts).at(0);
                    A(i, 2) = Psi3Der(particles, i + 3 * parts).at(0);
                }
                else if(coord == 2)
                {
                    A(i, 0) = Psi1Der(particles, i + 3 * parts).at(1);
                    A(i, 1) = Psi2Der(particles, i + 3 * parts).at(1);
                    A(i, 2) = Psi3Der(particles, i + 3 * parts).at(1);
                }
            }
        }
        else
        {
            A(i, 0) = Psi1(particles, i + 3 * parts);
            A(i, 1) = Psi2(particles, i + 3 * parts);
            A(i, 2) = Psi3(particles, i + 3 * parts);
        }
    }

    double det = A.determinant();

    return det;
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

std::vector<double> FermionsJastrow::Psi1Der(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
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

std::vector<double> FermionsJastrow::Psi2Der(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
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

std::vector<double> FermionsJastrow::Psi3Der(std::vector<std::unique_ptr<class Particle>>& particles, double part_inx)
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


std::vector<double> FermionsJastrow::gradpsi(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    std::vector<double> graddet(2);
    double denom = SD(particles, 4, 0, 0, 0);
    if(part_idx < 3)
    {
        graddet.at(0) = SD(particles, part_idx, 0, 1, 1)/denom;
        graddet.at(1) = SD(particles, part_idx, 0, 1, 2)/denom;
    }
    else
    {
        graddet.at(0) = SD(particles, part_idx - 3, 1, 1, 1)/denom;
        graddet.at(1) = SD(particles, part_idx - 3, 1, 1, 2)/denom;
    }

    return graddet;
}


std::vector<double> FermionsJastrow::gradlogjast(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    std::vector<double> gradjast(2);
    int N = 6;

    Particle& particle_j = *(particles.at(part_idx));
    double x_j = particle_j.getPosition().at(0);
    double y_j = particle_j.getPosition().at(1);

    double rij;
    double sum_1 = 0;
    double sum_2 = 0;
    for(int i = 0; i < 6; i++)
    {
        if(i == part_idx) continue;

        Particle particle_i = *(particles.at(i));
        double x_i = particle_i.getPosition().at(0);
        double y_i = particle_i.getPosition().at(1);

        rij = r_ij(particles, i, part_idx);
        
        int idx;
        if (i < part_idx)
        {
            idx = i * (2 * N - i - 1) / 2 + (part_idx - i - 1);
        }
        else
        {
            idx = part_idx * (2 * N - part_idx - 1) / 2 + (i - part_idx - 1);
        }

        double beta_ij = m_parameters.at(idx);

        sum_1 += beta_ij * (x_i - x_j) / rij;
        sum_2 += beta_ij * (y_i - y_j) / rij;
    }

    gradjast.at(0) = sum_1;
    gradjast.at(1) = sum_2;

    return gradjast;
}

double FermionsJastrow::gradpsijast(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double grad = 0;
    for(int i = 0; i < 6; i++)
    {
        std::vector<double> grad_jast = gradlogjast(particles, i);
        std::vector<double> grad_psi = gradpsi(particles, i);

        grad += grad_jast.at(0) * grad_psi.at(0) + grad_jast.at(1) * grad_psi.at(1);
    }

    return 2 * grad;
}

double FermionsJastrow::lapllogjast(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double laplacian = 0.0;
    double grad_squared = 0.0;
    int N = particles.size();

    for (int k = 0; k < N; ++k)
    {
        Particle& pk = *(particles[k]);
        double xk = pk.getPosition().at(0);
        double yk = pk.getPosition().at(1);

        double grad_k_x = 0.0;
        double grad_k_y = 0.0;

        for (int i = 0; i < N; ++i)
        {
            if (i == k) continue;

            Particle& pi = *(particles[i]);
            double xi = pi.getPosition().at(0);
            double yi = pi.getPosition().at(1);

            double dx = xk - xi;
            double dy = yk - yi;
            double r2 = dx * dx + dy * dy;
            double r = std::sqrt(r2);

            // Get index for beta_ik (symmetric, stored upper triangle)
            int idx = (i < k)
                ? i * (2 * N - i - 1) / 2 + (k - i - 1)
                : k * (2 * N - k - 1) / 2 + (i - k - 1);

            double beta_ik = m_parameters.at(idx);

            // Add to Laplacian ∇² log J
            laplacian += beta_ik * (1.0 / r - r2 / (r * r * r));

            grad_k_x += beta_ik * dx / r;
            grad_k_y += beta_ik * dy / r;
        }

        grad_squared += grad_k_x * grad_k_x + grad_k_y * grad_k_y;
    }

    return 2.0 * laplacian + grad_squared;
}

double FermionsJastrow::jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int N = particles.size();

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

double FermionsJastrow::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double det_up = SD(particles, 4, 0, 2, 0);
    double det_down = SD(particles, 4, 1, 2, 0);
    double J = jastrow(particles);

    //std:: cout << det_up * det_down << std::endl;

    return det_up * det_down * J;
}

double FermionsJastrow::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double laplacian1 = 0;
    double laplacian2 = 0;
    for (int i = 0; i < 3 ; i++)
    {
        laplacian1 += SD(particles, i, 0, 2, 0);
        laplacian2 += SD(particles, i, 1, 2, 0);
    }
    laplacian1 /=  SD(particles, 4, 0, 2, 0);
    laplacian2 /=  SD(particles, 4, 1, 2, 0);

    double laplacianpsi = (laplacian1 + laplacian2);
    double grad_psijast = gradpsijast(particles);
    double laplacian_jast = lapllogjast(particles);

    //std::cout << laplacianpsi << std::endl;
    //std::cout << grad_psijast << std::endl;
    //std::cout << laplacian_jast << std::endl;


    return laplacianpsi + grad_psijast + laplacian_jast;
}