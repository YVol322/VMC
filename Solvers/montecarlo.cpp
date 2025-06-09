#include "montecarlo.h"     // Include "montecarlo" header file with declarations.


// Constructor of the MonteCarlo class. It initializes the Monte Carlo simulation with a random number generator.
//
// Input:   std::unique_ptr<class Random> rng - a unique pointer to a random number generator instance.
//
// Output:  void - no return value.
MonteCarlo::MonteCarlo(std::unique_ptr<class Random> rng)
{
    m_rng = std::move(rng);
}