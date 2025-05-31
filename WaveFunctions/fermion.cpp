#include "fermion.h"

Fermion::Fermion(double alpha, std::vector<double> beta, int mode, int n_particles)
{

    int n_betas = beta.size();
    m_numberOfParameters = n_betas + 1;

    m_parameters.reserve(m_numberOfParameters);
    m_parameters.insert(m_parameters.end(), beta.begin(), beta.end());
    m_parameters.push_back(alpha);

    m_particles = n_particles;
    m_mode = mode;
}


double Fermion::Jastrow(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double sum = 0.0;
    int n_particles = m_particles;

    for (int i = 0; i < n_particles - 1; i++)
    {
        for (int j = i + 1; j < n_particles; j++)
        {
            double rij = r_ij(particles, i, j);

            int idx = BetaIndex(i,j);
            double beta_ij = m_parameters[idx];

            sum += beta_ij * rij;
        }
    }

    return exp(sum);
}


double Fermion::PadeJastrow(std::vector<std::unique_ptr<class Particle>>& particles)
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


double Fermion::Psi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    double r2 = r_squared(particles, part_idx);

    return exp(-alpha * r2);
}


double Fermion::Psi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double x = particles[part_idx] -> getPosition()[0];

    return x * Psi1(particles, part_idx);
}

double Fermion::Psi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double y = particles[part_idx] -> getPosition()[1];

    return y * Psi1(particles, part_idx);
}


double Fermion::Psi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];


    return x * y * Psi1(particles, part_idx);
}


double Fermion::Psi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double x = particles[part_idx] -> getPosition()[0];

    return (x * x - 1) * Psi1(particles, part_idx);
}


double Fermion::Psi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double y = particles[part_idx] -> getPosition()[1];

    return (y * y - 1) * Psi1(particles, part_idx);
}

std::vector<double> Fermion::GradiPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi1 = {-2 * alpha * x * psi1, -2 * alpha * y * psi1};

    return grad_psi1;
}

std::vector<double> Fermion::GradiPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi2 = {(1 - 2 * alpha * x * x) * psi1, -2 * alpha * y * x * psi1};

    return grad_psi2;
}

std::vector<double> Fermion::GradiPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi3 = {-2 * alpha * x * y * psi1, (1 - 2 * alpha * y * y) * psi1};

    return grad_psi3;
}


std::vector<double> Fermion::GradiPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi4 = {y * psi1 * (1 - 2 * alpha * x * x), x * psi1 * (1 - 2 * alpha * y * y)};

    return grad_psi4;
}


std::vector<double> Fermion::GradiPsi5(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi5 = {2 * x * psi1 * (1 - alpha * x * x + alpha), -2 * alpha * y * (x * x - 1) * psi1};

    return grad_psi5;
}


std::vector<double> Fermion::GradiPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    std::vector<double> grad_psi6 = {-2 * alpha * x * (y * y - 1) * psi1, 2 * y * psi1 * (1 - alpha * y * y + alpha)};

    return grad_psi6;
}


double Fermion::LapliPsi1(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-4 * alpha + 4 * alpha * alpha * (x * x + y * y) ) * Psi1(particles, part_idx);
}

double Fermion::LapliPsi2(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-8 * alpha * x + 4 * alpha * alpha * x * (x * x + y * y) ) * Psi1(particles, part_idx);
}

double Fermion::LapliPsi3(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return (-8 * alpha * y + 4 * alpha * alpha * y * (x * x + y * y) ) * Psi1(particles, part_idx);
}


double Fermion::LapliPsi4(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();
    
    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];

    return x * y * Psi1(particles, part_idx) * (4 * alpha * alpha * (x * x + y * y) - 12 * alpha);
}


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
    (-4 * alpha * x2 + alpha * common * (2 * alpha * x2 - 1) + alpha * common * (2 * alpha * y2 - 1) + 1);

    return laplacian_psi5;
}


double Fermion::LapliPsi6(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    double alpha = m_parameters.back();

    double x = particles[part_idx] -> getPosition()[0];
    double y = particles[part_idx] -> getPosition()[1];
    double psi1 = Psi1(particles, part_idx);

    double x2 = x * x;
    double y2 = y * y;
    double common = y2 - 1;

    double laplacian_psi6 = 2 * psi1 *
    (-4 * alpha * y2 + alpha * common * (2 * alpha * x2 - 1) + alpha * common * (2 * alpha * y2 - 1) + 1);

    return laplacian_psi6;
}


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


std::vector<double> Fermion::GradiJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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


std::vector<double> Fermion::GradiPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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


double Fermion::LapliJOverJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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


double Fermion::LapliPJOverPJ(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
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


double Fermion::LaplPsi1OverPsi1(std::vector<std::unique_ptr<class Particle>>& particles)
{

    double det_up = SD(particles, 0, 0, 0, 0);
    double det_down = SD(particles, 1, 0, 0, 0);

    double sum = 0.0;
    for(int i = 0; i < m_particles / 2; i++)
    {
        sum += SD(particles, 0, i, 2, 0) / det_up;
        sum += SD(particles, 1, i, 2, 0) / det_down;
    }

    return sum;
}


double Fermion::GradiPsi1GradiJOverPsi(std::vector<std::unique_ptr<class Particle>>& particles, double part_idx)
{
    int n_dimensions = particles[0] -> getNumberOfDimensions();

    std::vector<double> grad_J = (m_mode == 0) ? GradiJOverJ(particles, part_idx) : GradiPJOverPJ(particles, part_idx);

    double det_up = SD(particles, 0, 0, 0, 0);
    double det_down = SD(particles, 1, 0, 0, 0);

    double sum = 0.0;
    if (part_idx < m_particles/2)
    {
        for (int i = 0; i < n_dimensions; i++)
        {
            sum += SD(particles, 0, part_idx, 1, i) / det_up * grad_J[i];
        }
    }
    else
    {
        int idx = part_idx - m_particles/2;
        for (int i = 0; i < n_dimensions; i++)
        {
            sum += SD(particles, 1, idx, 1, i) / det_down * grad_J[i];
        }
    }

    return 2.0 * sum;
}


double Fermion::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double J = (m_mode == 0) ? Jastrow(particles) : PadeJastrow(particles);

    return SD(particles, 0, 0, 0, 0) * SD(particles, 1, 0, 0, 0) * J;
}

double Fermion::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double LaplJ = 0;
    double GradGrad = 0;

    if(m_mode == 0)
    {
        for(int i = 0; i < m_particles; i++)
        {
            LaplJ += LapliJOverJ(particles, i);
            GradGrad += GradiPsi1GradiJOverPsi(particles, i);
        }
    }
    else
    {
        for(int i = 0; i < m_particles; i++)
        {
            LaplJ += LapliPJOverPJ(particles, i);
            GradGrad += GradiPsi1GradiJOverPsi(particles, i);
        }
    }

    return LaplPsi1OverPsi1(particles) + LaplJ + GradGrad;
}


std::vector<double> Fermion::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int i)
{
    int half = m_particles/2;
    int spinBlock = (i < half ? 0 : 1);
    int localIdx  = (i < half ? i : i - half);

    double dlnD_dx = SD(particles, spinBlock, localIdx, 1, 0) / SD(particles, spinBlock, 0, 0, 0);
    double dlnD_dy = SD(particles, spinBlock, localIdx, 1, 1) / SD(particles, spinBlock, 0, 0, 0);

    auto dlnJ = (m_mode == 0) ? GradiJOverJ(particles, i) : GradiPJOverPJ(particles, i);

    std::vector<double> F = { 2*(dlnD_dx + dlnJ[0]), 2*(dlnD_dy + dlnJ[1]) };
    
    return F;
}