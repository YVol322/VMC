#include "bosonnumerical.h"


BosonNumerical::BosonNumerical(double alpha, int n_particles, double omega)
{
    assert(alpha >= 0);

    m_numberOfParameters = 1;

    m_parameters.reserve(m_numberOfParameters);
    m_parameters.push_back(alpha);

    m_particles = n_particles;

    m_omega = omega;
    m_sqrt_om = sqrt(omega);
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


var BosonNumerical::PsiT(VectorXvar& x)
{
    int n_dimensions = x.size() / m_particles;

    var sum = 0;
    for(int i = 0; i < m_particles; i++)
    {
        for(int j = 0; j < n_dimensions; j++)
        {
            var pos = x(i * n_dimensions + j);
            sum += pos * pos;
        }
    }
    var alpha = m_parameters.back();

    return exp(-alpha * m_omega * sum);
}


double BosonNumerical::LaplPsiTOverPsiT(VectorXvar& x)
{
    var psiT = PsiT(x);

    Eigen::VectorXd g;
    Eigen::MatrixXd H = hessian(psiT, x, g);
    
    Eigen::VectorXd LaplVector= H.diagonal();

    double sum = val(LaplVector.sum());

    return sum / val(psiT);
}



double BosonNumerical::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    VectorXvar x = fill_x(particles);

    return val(PsiT(x));
}


double BosonNumerical::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    VectorXvar x = fill_x(particles);

    return LaplPsiTOverPsiT(x);
}


std::vector<double> BosonNumerical::quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int n)
{
    VectorXvar x = fill_x(particles);
    var psiT = PsiT(x);

    int d = x.size() / m_particles;
    std::vector<double> F(d);
    
    for(int i = 0; i < d; i++)
    {
        auto [gradx] = derivatives(psiT, wrt(x(d*n + i)));
        F[i] = 2 * val(gradx) / (val(psiT));
    }
    
    return F;
}