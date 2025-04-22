#include <iostream>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

using namespace autodiff;

real sp(const ArrayXreal& x, int idx, real alpha)
{
    return exp(-alpha * (x(idx) * x(idx) + x(idx + 1) * x(idx + 1)));
}

real jastrow(const ArrayXreal& x, real beta)
{
    return exp(-beta * sqrt( (x(0) - x(2)) * (x(0) - x(2)) + (x(1) - x(3)) * (x(1) - x(3)) ));
}

// The scalar function for which the gradient is needed
real f(const ArrayXreal& x, real beta, real alpha)
{
    real psi = 1;
    for(int i = 0; i < 2; i++)
    {
        psi *= sp(x, i, alpha);
    }
    psi *= jastrow(x, beta);

    return psi;
}

int main()
{
    using Eigen::VectorXd;

    real alpha = 0.5;
    real beta = 0;

    ArrayXreal x(4);                            // the input array x with 4 variables
    x << 0.238754, 0.473744, 0.132144, 0.439087;                         // x = [1, 2, 3, 4]

    real u;                                     // the output scalar u = f(x) evaluated together with gradient below

    // Compute the gradient (du/dx)
    VectorXd g = gradient(f, wrt(x), at(x, beta, alpha), u); // first gradient, du/dx

    // Now compute the gradient of the gradient (i.e., second derivatives)
    // This can be thought of as the Laplacian, summing the second derivatives.
    VectorXd laplacian(g.size());
    for (int i = 0; i < g.size(); ++i)
    {
        // Compute the second derivative by differentiating the gradient with respect to x(i)
        laplacian(i) = derivative(f, wrt(x(i)), at(x, beta, alpha), u); // second derivative, d^2u/dx_i^2
    }

    // Sum the second derivatives to get the Laplacian
    real lap = laplacian.sum();

    std::cout << "u = " << u << std::endl;         // print the evaluated output u
    std::cout << "Gradient (du/dx) = \n" << g << std::endl;    // print the gradient vector
    std::cout << "Laplacian (sum of second derivatives) = " << lap << std::endl;  // print the Laplacian

    return 0;
}

// -2.66344 -2.66344 -2.66344 -2.66344