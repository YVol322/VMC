#include "boson.h"


Boson::Boson(double alpha, std::vector<double> beta, int mode, int n_particles)
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


double Boson::Phi(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    int n_particles = m_particles;
    double argument = 0;

    for(int i = 0; i < n_particles; i++)
    {
        argument += r_squared(particles, i);
    }

    return exp(-alpha * argument);
}


int Boson::BetaIndex(int i, int j)
{
    int n_particles = m_particles;

    int ii = std::min(i, j);
    int jj = std::max(i, j);

    return ii * (2 * n_particles - ii - 1) / 2 + (jj - ii - 1);
}


double Boson::a_ij(int i, int j)
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


double Boson::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;

    for (int i = 0; i < n_particles - 1; ++i)
    {
        for (int j = i + 1; j < n_particles; ++j)
        {
            double rij = r_ij(particles, i, j);

            int idx = BetaIndex(i,j);
            double beta_ij = m_parameters[idx];

            sum += beta_ij * rij;
        }
    }

    return exp(sum);
}


double Boson::PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;
    double beta = m_parameters[0];

    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {
            double aij = a_ij(i, j);
            double rij = r_ij(particles, i, j);

            sum += aij * rij / (1 + beta * rij);
        }
    }

    return exp(sum);
}


std::vector<double> Boson::GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    Particle& particle_i = *(particles[part_idx]);

    std::vector<double> grad(n_dimensions, 0.0);

    for(int j = 0; j < n_particles; j++)
    {
        if(j == part_idx) continue;

        Particle& particle_j = *(particles[j]);
        double ril = r_ij(particles, part_idx, j);

        int idx = BetaIndex(part_idx, j);
        double beta_il = m_parameters[idx];

        for(int d = 0; d < n_dimensions; d++)
        {
            double coord_i = particle_i.getPosition()[d];
            double coord_j = particle_j.getPosition()[d];

            grad[d] += beta_il * (coord_i - coord_j) / ril;
        }
    }

    return grad;
}


std::vector<double> Boson::GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    int n_particles = particles.size();
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    double beta = m_parameters.at(0);

    Particle particle_i = *(particles[part_idx]);

    std::vector<double> grad(n_dimensions, 0.0);

    for(int j = 0; j < n_particles; j++)
    {
        if(j == part_idx) continue;

        Particle particle_j = *(particles[j]);
        double ril = r_ij(particles, part_idx, j);
        double aij = a_ij(part_idx, j);

        for(int d = 0; d < n_dimensions; d++)
        {
            double coord_i = particle_i.getPosition()[d];
            double coord_j = particle_j.getPosition()[d];

            grad[d] += aij * (coord_i - coord_j) / (ril * (1 + beta * ril) * (1 + beta * ril));
        }
    }

    return grad;
}


double Boson::LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    std::vector<double> grad_i = GradiJOverJ(particles, part_idx);
    double sum1 = 0.0;
    for (int d = 0; d < n_dimensions; d++)
    {
        sum1 += grad_i[d] * grad_i[d];
    }

    double sum2 = 0.0;
    for (int j = 0; j < n_particles; j++)
    {
        if (j == part_idx) continue;

        int idx = BetaIndex(part_idx, j);
        double beta_il = m_parameters[idx];

        double ril = r_ij(particles, part_idx, j);

        sum2 += beta_il / ril;
    }

    return sum1 + sum2;
}


double Boson::LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
{
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    std::vector<double> grad_i = GradiPJOverPJ(particles, part_idx);
    double sum1 = 0.0;
    for (int d = 0; d < n_dimensions; d++)
    {
        sum1 += grad_i[d] * grad_i[d];
    }

    double sum2 = 0.0;
    double beta = m_parameters[0];

    for (int j = 0; j < n_particles; j++)
    {
        if (j == part_idx) continue;

        double aij = a_ij(part_idx, j);
        double ril = r_ij(particles, part_idx, j);
        double t = 1 + beta * ril;

        sum2 += aij * (1/(ril * t * t) - 2 * beta/(t * t * t));
    }

    return sum1 + sum2;
}


double Boson::LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double alpha = m_parameters.back();
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

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
    int n_particles = m_particles;
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    double sum1 = 0;
    for(int i = 0; i < n_particles; i++)
    {
        Particle& particle_i = *(particles[i]);
        for(int j = 0; j < n_particles; j++)
        {
            if (j == i) continue;

            Particle& particle_j = *(particles[j]);

            double rij = r_ij(particles, i, j);

            double sum2 = 0;
            for(int d = 0; d < n_dimensions; d++)
            {
                double coord_i = particle_i.getPosition()[d];
                double coord_j = particle_j.getPosition()[d];
                sum2 += coord_i * (coord_i - coord_j);
            }
            if(m_mode == 0) 
            {
                int idx = BetaIndex(i, j);
                double beta_ij = m_parameters.at(idx);

                sum1 += sum2 * beta_ij / rij;
            }
            else 
            {
                double beta = m_parameters.at(0);
                double aij = a_ij(i, j);
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
    if(m_mode == 0) return Phi(particles) * Jastrow(particles);
    else return Phi(particles) * PadeJastrow(particles);
}


double Boson::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
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


double Boson::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
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


double Boson::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
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