#include "particle.h"   // Include "particle" header file with declarations.



// Constructor of the Particle class. It takes a const std::vector<double>& position as its argument,
// sets the private variable m_position equal to it, and saves the size of this vector to the private variable
// m_numberOfDimensions.
Particle::Particle(const std::vector<double>& position)
{
    m_numberOfDimensions = position.size();
    m_position = position;
}



// Function that adjusts the position of the particle.
//
// Input:   double change - the adjustment value;
//          unsigned int dimension - specifies which coordinate is changed.
//
// Output: void - no return value.
void Particle::adjustPosition(double change, unsigned int dimension)
{
    m_position.at(dimension) += change;
}