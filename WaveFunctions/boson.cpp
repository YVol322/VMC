#include <cmath>
#include <iostream>

#include "boson.h"

Boson::Boson(double alpha, std::vector<double> beta, int mode, int n_particles)
{
    assert(alpha >= 0);
    int n_betas = beta.size();
    m_numberOfParameters = n_betas + 1;
    m_parameters.reserve(m_numberOfParameters);

    for(int i = 0; i < n_betas; i++)
    {
        m_parameters.push_back(beta.at(i));
    }

    m_particles = n_particles;
    m_mode = mode;

    m_parameters.push_back(alpha);
}


double Boson::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;

    for (int i = 0; i < n_particles - 1; ++i)
    {
        for (int j = i + 1; j < n_particles; ++j)
        {
            double rij = r_ij(particles, i, j);

            int idx = i * (2 * n_particles - i - 1) / 2 + (j - i - 1);
            double beta_ij = m_parameters.at(idx);

            sum += beta_ij * rij;
        }
    }

    return std::exp(sum);
}


double Boson::PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;
    double aij = 0;

    for (int i = 0; i < n_particles - 1; ++i)
    {
        for (int j = i + 1; j < n_particles; ++j)
        {
            if(i < n_particles/2 && j >= n_particles/2)
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


std::vector<double> Boson::GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0]->getNumberOfDimensions();

    Particle& particle_i = *(particles.at(part_inx));


    std::vector<double> grad(n_dimensions, 0.0);

    for(int l = 0; l < n_particles; l++)
    {
        if(l == part_inx) continue;

        Particle& particle_l = *(particles.at(l));
        double ril = r_ij(particles, part_inx, l);

        int i = std::min(part_inx, l);
        int j = std::max(part_inx, l);
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


std::vector<double> Boson::GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx)
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


double Boson::LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx)
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

        int i = std::min(part_inx, l);
        int j = std::max(part_inx, l);
        int idx = i * (2 * n_particles - i - 1) / 2 + (j - i - 1);
        double beta_il = m_parameters.at(idx);

        double ril = r_ij(particles, part_inx, l);

        sum2 += beta_il / ril;
    }

    return sum1 + sum2;
}


double Boson::LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_inx)
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


double Boson::LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    double n_particles = m_particles;
    double n_dimensions = particles[0] -> getNumberOfDimensions();

    double r2 = 0;
    for(int i = 0; i < n_particles; i++)
    {
        r2 += r_squared(particles, i);
    }

    return 2 * alpha * (2 * alpha * r2 - n_particles * n_dimensions);
}


double Boson::GradPsi1GradJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    double n_particles = m_particles;
    double n_dimensions = particles[0] -> getNumberOfDimensions();

    double sum1 = 0;
    for(int i = 0; i < n_particles; i++)
    {
        Particle& particle_i = *(particles.at(i));
        for(int j = 0; j < n_particles; j++)
        {
            if (j == i) continue;

            Particle& particle_j = *(particles.at(j));

            double rij = r_ij(particles, i, j);

            double sum2 = 0;
            for(int k = 0; k < n_dimensions; k++)
            {
                double coord_i = particle_i.getPosition().at(k);
                double coord_j = particle_j.getPosition().at(k);
                sum2 += coord_i * (coord_i - coord_j);
            }
            if(m_mode == 0) 
            {
                int ii = std::min(i, j);
                int jj = std::max(i, j);
                int idx = ii * (2 * n_particles - ii - 1) / 2 + (jj - ii - 1);
                double beta_ij = m_parameters.at(idx);
                sum1 += sum2 * beta_ij / rij;
            }
            else 
            {
                double beta = m_parameters.at(0);
                double aij = (i < n_particles / 2 && j >= n_particles / 2) || (i >= n_particles / 2 && j < n_particles / 2)
                             ? 1.0 : 1.0 / 3.0;
                double t = 1 + beta * rij;
                sum1 += sum2 * aij / (t * t * rij);
            }
        }
    }

    sum1 *= -4 * alpha;

    return sum1;
}


double Boson::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    double alpha = m_parameters.back();
    int n_particles = m_particles;
    double argument = 0;
    double wavefunction;

    for(int i = 0; i < n_particles; i++)
    {
        argument += r_squared(particles, i);
    }

    argument *= -alpha;
    wavefunction = exp(argument);

    if(m_mode == 0 ) return wavefunction * Jastrow(particles);
    else return wavefunction * PadeJastrow(particles);
}

double Boson::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    int n_particles = m_particles;

    double lapl_Psi1 = LaplPsi1OverPsi1(particles);
    double gradgrad = GradPsi1GradJOverPsi(particles);
    double lapl_J = 0;
    if(m_mode == 0)
    {
        for(int i = 0; i < n_particles; i++)
        {
            lapl_J += LapliJOverJ(particles, i);
        }
    }
    else
    {
        for(int i = 0; i < n_particles; i++)
        {
            lapl_J += LapliPJOverPJ(particles, i);
        }
    }

    return lapl_Psi1 + gradgrad + lapl_J;
}


double Boson::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
{
    Particle& particle_i = *(particles.at(part_index));
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


double Boson::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
{
    Particle& pi = *(particles.at(i));
    Particle& pj = *(particles.at(j));

    int n_dimensions = pi.getNumberOfDimensions();

    double rij = 0;
    for(int d = 0; d < n_dimensions; d++)
    {
        rij += (pi.getPosition().at(d) - pj.getPosition().at(d)) * (pi.getPosition().at(d) - pj.getPosition().at(d));
    }
    
    return sqrt(rij);
}