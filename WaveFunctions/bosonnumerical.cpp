#include <cmath>
#include <iostream>

#include "bosonnumerical.h"

BosonNumerical::BosonNumerical(double alpha, std::vector<double> beta, int mode, int n_particles)
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


VectorXvar BosonNumerical::fill_x(std::vector<std::unique_ptr<class Particle>>& particles)
{
    int n_particles = m_particles;
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


var BosonNumerical::Jastrow(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;

    var sum = 0.0;
    for (int i = 0; i < n_particles - 1; ++i)
    {
        for (int j = i + 1; j < n_particles; ++j)
        {
            var rij2 = 0.0;
            for (int d = 0; d < n_dimensions; d++)
            {
                var delta = x(i * n_dimensions + d) - x(j * n_dimensions + d);
                rij2 += delta * delta;
            }
            var rij = sqrt(rij2);

            int idx = i * (2 * n_particles - i - 1) / 2 + (j - i - 1);
            var beta_ij = m_parameters.at(idx);

            sum += beta_ij * rij;
        }
    }

    return exp(sum);
}


var BosonNumerical::PadeJastrow(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;
    var sum = 0.0;

    for (int i = 0; i < n_particles - 1; ++i)
    {
        for (int j = i + 1; j < n_particles; ++j)
        {
            var rij2 = 0.0;
            for (int d = 0; d < n_dimensions; d++)
            {
                var delta = x(i * n_dimensions + d) - x(j * n_dimensions + d);
                rij2 += delta * delta;
            }
            var rij = sqrt(rij2);

            var beta = m_parameters.at(0);

            var aij;
            if(i < n_particles/2 && j >= n_particles/2)
            {
                aij = 1.0;
            }
            else
            {
                aij = 1.0/3.0;
            }

            sum += aij * rij / (1 + beta * rij);
        }
    }

    return exp(sum);
}


var BosonNumerical::Phi(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;
    var sum = 0;
    for(int i = 0; i < n_particles; i++)
    {
        for(int j = 0; j < n_dimensions; j++)
        {
            sum += x(i * n_dimensions + j) * x(i * n_dimensions + j);
        }
    }
    var alpha = m_parameters.back();

    sum *= -alpha;

    return exp(sum);
}


VectorXvar BosonNumerical::dJdxi(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;

    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);
    VectorXvar dJ = gradient(J, x);

    return dJ;
}


VectorXvar BosonNumerical::dPhidxi(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;

    var phi = Phi(x);
    VectorXvar dPhi = gradient(phi, x);

    return dPhi;
}


VectorXvar BosonNumerical::d2Phidxi2(VectorXvar& x)
{
    var u = Phi(x);

    Eigen::VectorXd g;
    Eigen::MatrixXd H = hessian(u, x, g);
    
    return H.diagonal();
}


VectorXvar BosonNumerical::d2Jdxi2(VectorXvar& x)
{
    var u = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);

    Eigen::VectorXd g;
    Eigen::MatrixXd H = hessian(u, x, g);
    
    return H.diagonal();
}


double BosonNumerical::LaplJOverJ(VectorXvar& x)
{
    int n = x.size();
    VectorXvar d2J(n);
    d2J = d2Jdxi2(x);

    double sum = val(d2J.sum());

    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);

    return sum / val(J);
}


double BosonNumerical::LaplPhiOverPhi(VectorXvar& x)
{
    int n = x.size();
    VectorXvar d2Phi(n);
    d2Phi = d2Phidxi2(x);

    double sum = val(d2Phi.sum());
    double phi = val(Phi(x));

    return sum / phi;
}


double BosonNumerical::GradGrad(VectorXvar& x)
{
    int n = x.size();
    VectorXvar dPhi(n), dJ(n);
    dPhi = dPhidxi(x);
    dJ = dJdxi(x);

    double phi = val(Phi(x));
    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);


    return 2.0 * val(dJ.dot(dPhi)) / (val(J) * phi);
}


double BosonNumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    VectorXvar x = fill_x(particles);

    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);
    double phi = val(Phi(x));

    return val(J) * phi;
}

double BosonNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);
    
    double LaplJ = LaplJOverJ(x);
    double LaplPhi = LaplPhiOverPhi(x);
    double gradgrad = GradGrad(x);

    return LaplJ + gradgrad + LaplPhi;
}


double BosonNumerical::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_index)
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


double BosonNumerical::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
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