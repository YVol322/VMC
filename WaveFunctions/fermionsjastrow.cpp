#include <cmath>
//#include <armadillo>
#include <Eigen/Dense>
#include <iostream>


#include "fermionsjastrow.h"

FermionsJastrow::FermionsJastrow(double alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(m_numberOfParameters);
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

double FermionsJastrow::SD(std::vector<std::unique_ptr<class Particle>>& particles, double row_indx, double parts)
{
    Eigen::MatrixXd A(3, 3);

    for(int i = 0; i < 3; i++)
    {
        if(i == row_indx)
        {
            A(i, 0) = Psi1DoubleDer(particles, i + 3 * parts);
            A(i, 1) = Psi2DoubleDer(particles, i + 3 * parts);
            A(i, 2) = Psi3DoubleDer(particles, i + 3 * parts);
        }
        else
        {
            A(i, 0) = Psi1(particles, i + 3 * parts);
            A(i, 1) = Psi2(particles, i + 3 * parts);
            A(i, 2) = Psi3(particles, i + 3 * parts);
        }

        //std::cout << A(i, 0) << std::endl;
        //std::cout << A(i, 1) << std::endl;
        //std::cout << A(i, 2) << std::endl;
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

double FermionsJastrow::evaluate(std::vector<std::unique_ptr<class Particle>>& particles)
{   
    double det_up = SD(particles, 4, 0);
    double det_down = SD(particles, 4, 1);

    //std:: cout << det_up * det_down << std::endl;

    return det_up * det_down;
}

double FermionsJastrow::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles)
{
    double laplacian1 = 0;
    double laplacian2 = 0;
    for (int i = 0; i < 3 ; i++)
    {
        laplacian1 += SD(particles, i, 0);
        laplacian2 += SD(particles, i, 1);
    }
    laplacian1 /=  SD(particles, 4, 0);
    laplacian2 /=  SD(particles, 4, 1);

    double laplacian = (laplacian1 + laplacian2);

    //std::cout << laplacian << std::endl;


    return laplacian;
}