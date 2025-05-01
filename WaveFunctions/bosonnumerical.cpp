#include "bosonnumerical.h"


BosonNumerical::BosonNumerical(double alpha, std::vector<double> beta, int mode, int n_particles)
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


VectorXvar BosonNumerical::fill_x(std::vector<std::unique_ptr<class Particle>>& particles)
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


int BosonNumerical::BetaIndex(int i, int j)
{
    int n_particles = m_particles;

    int ii = std::min(i, j);
    int jj = std::max(i, j);

    return ii * (2 * n_particles - ii - 1) / 2 + (jj - ii - 1);
}


var BosonNumerical::a_ij(int i, int j)
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


var BosonNumerical::Jastrow(VectorXvar& x)
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


var BosonNumerical::PadeJastrow(VectorXvar& x)
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


var BosonNumerical::Phi(VectorXvar& x)
{
    int n_particles = m_particles;
    int n_dimensions = x.size() / n_particles;

    var sum = 0;
    for(int i = 0; i < n_particles; i++)
    {
        for(int j = 0; j < n_dimensions; j++)
        {
            var pos = x(i * n_dimensions + j);
            sum += pos * pos;
        }
    }
    var alpha = m_parameters.back();

    return exp(-alpha * sum);
}


VectorXvar BosonNumerical::dJdxi(VectorXvar& x)
{
    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);
    VectorXvar dJ = gradient(J, x);

    return dJ;
}


VectorXvar BosonNumerical::dPhidxi(VectorXvar& x)
{
    var phi = Phi(x);
    VectorXvar dPhi = gradient(phi, x);

    return dPhi;
}


VectorXvar BosonNumerical::d2Phidxi2(VectorXvar& x)
{
    var phi = Phi(x);

    Eigen::VectorXd g;
    Eigen::MatrixXd H = hessian(phi, x, g);
    
    return H.diagonal();
}


VectorXvar BosonNumerical::d2Jdxi2(VectorXvar& x)
{
    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);

    Eigen::VectorXd g;
    Eigen::MatrixXd H = hessian(J, x, g);
    
    return H.diagonal();
}


double BosonNumerical::LaplJOverJ(VectorXvar& x)
{
    VectorXvar d2J = d2Jdxi2(x);

    double sum = val(d2J.sum());

    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);

    return sum / val(J);
}


double BosonNumerical::LaplPhiOverPhi(VectorXvar& x)
{
    VectorXvar d2Phi = d2Phidxi2(x);

    double sum = val(d2Phi.sum());

    return sum / val(Phi(x));
}


double BosonNumerical::GradGrad(VectorXvar& x)
{
    VectorXvar dPhi = dPhidxi(x);
    VectorXvar dJ = dJdxi(x);

    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);

    return 2.0 * val(dJ.dot(dPhi)) / (val(J) * val(Phi(x)));
}


double BosonNumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    VectorXvar x = fill_x(particles);

    var J = (m_mode == 0) ? Jastrow(x) : PadeJastrow(x);

    return val(J) * val(Phi(x));
}


double BosonNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    return LaplJOverJ(x) + LaplPhiOverPhi(x) + GradGrad(x);
}


double BosonNumerical::r_squared(std::vector<std::unique_ptr<class Particle>>& particles, int part_idx)
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


double BosonNumerical::r_ij(std::vector<std::unique_ptr<class Particle>>& particles, int i, int j)
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