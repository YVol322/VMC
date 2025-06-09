# Yevhenii Volkov's Variational Monte Carlo Solver for many-body systems.


## Table of Contents

- [Description](#description)
- [Features](#features)
- [Installation](#installation)
- [Compilation and Execution Steps](#compilation-and-execution-steps)
- [Configuration](#configuration)
- [Examples](#examples)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [Contributing](#contributing)


## Description

This directory contains a flexible, object-oriented C++ implementation of a Variational Monte Carlo (VMC) solver designed to estimate the ground-state energy of selected many-body systems. The primary systems studied include:

- **Non-interacting bosons**: \(N\) non-interacting bosons in \(d\) dimensions confined in a harmonic oscillator (HO) potential, modeling a Bose-Einstein condensate.
- **Non-interacting electrons**: \(2\), \(6\), and \(12\) non-interacting electrons in \(2\) dimensions confined in a HO potential, modeling an idealized 2D quantum dot.
- **Interacting electrons**: \(2\), \(6\), and \(12\) interacting electrons in \(2\) dimensions confined in a HO potential with Coulomb interaction, modeling a more realistic 2D quantum dot.

### Hamiltonian Formulations

- **Non-interacting particles**:  
  \[
  \hat{H} = \sum_{i=1}^N \left[ -\frac{1}{2}\nabla_i^2 + \frac{1}{2} m\omega^2 r_i^2 \right]
  \]

- **Interacting electrons**:  
  \[
  \hat{H} = \sum_{i=1}^N \left[ -\frac{1}{2}\nabla_i^2 + \frac{1}{2} m\omega^2 r_i^2 \right] + \sum_{i=1}^N \sum_{j>i}^N \frac{1}{r_{ij}}
  \]

For the interacting electrons, the system includes two types of correlation factors:
- **Jastrow factor**:  
  \[
  J = \exp \left( \sum_{i=1}^N \sum_{j>i}^N \beta_{ij} r_{ij} \right)
  \]
  
- **Pade-Jastrow factor**:  
  \[
  P = \exp \left( \sum_{i=1}^N \sum_{j>i}^N \frac{a r_{ij}}{1 + \beta r_{ij}} \right)
  \]

### Wavefunctions and Variational Parameters

- **Non-interacting wavefunctions** are modified by a single variational parameter \(\alpha\).
- **Interacting Jastrow ansatz** involves \(p = \frac{N(N-1)}{2}\) variational parameters \(\beta_{ij}\), where \(N\) is the number of particles.
- **Pade-Jastrow ansatz** is modified with one variational parameter \(\beta\).

These parameters are optimized using Gradient Descent to minimize the variational energy, thus providing an estimate for the ground-state energy.

### Sampling Algorithms

The project implements two sampling algorithms:

1. **Metropolis Algorithm**: Uses uniform random numbers to adjust the positions of all particles simultaneously.
2. **Metropolis-Hastings Algorithm**: Uses Green's function as a proposal distribution, based on Fokker-Planck and Langevin equations. This method moves one particle at a time, accepts or rejects the step, and repeats for all particles.

### Automatic Differentiation and Parallelization

The VMC solver uses the **autodiff** C++ API to implement automatic differentiation. While this introduces some computational overhead, it ensures the correctness of the analytical expressions.

In addition, the main VMC programs are parallelized using **OpenMP** and **MPI** to reduce overhead and speed up the algorithms. The MPI implementation allows the code to be executed on supercomputers or cluster computers, enabling the simulation of larger systems.

### Flexibility and Extensibility

The structure of the code is highly flexible. Trial wavefunctions, sampling algorithms, and Hamiltonians can be easily modified or extended, making it straightforward to study other systems. New features can also be added as needed.

## Features

- **Flexible object-oriented C++ implementation** for solving variational Monte Carlo problems.
- **Various systems** studied, including non-interacting and interacting bosons/electrons in 2D.
- **Support for Jastrow and Pade-Jastrow correlation factors** for interacting particles.
- **Optimized using Gradient Descent** for variational energy minimization.
- **Parallelization with OpenMP and MPI** for efficient execution on multi-core and distributed systems.
- **Automatic differentiation** implemented using the autodiff C++ API.


## Installation

### Prerequisites

To compile and run the program, make sure you have the following:

- **C++17 compiler** (e.g., GCC, Clang)
- **OpenMP** for C++17
- **MPI** for C++17
- **Autodiff library** for automatic differentiation

### Installation Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/YVol322/VMC.git

## Compilation and Execution Steps

To compile and execute the program, follow these steps:

1. Copy the desired program from either the `MainFunctions` or `ResultsGeneration` directories to the main VMC directory.

2. To compile and run the **serial** or **OpenMP** version of the program, execute the following command:

   ```bash
   ./execute

3. To compile and run the MPI version of the program, execute:
   ```bash
    ./execute_MPI

This will compile and run the program depending on the selected execution method (serial, OpenMP, or MPI).

### Important Note
Due to the CMake setup, only one executable .cpp file can be present in the VMC directory at a time. If more than one executable .cpp file is found, the compilation will fail.


## Configuration

Before running the VMC solver, you may need to configure certain parameters. The following options can be set:

- **Number of particles**: Adjust the number of particles in the system. 
  - For fermions, only 2, 6, and 12 particles in 2D are supported.
  - For bosons, \(N\) particles in \(d\) dimensions can be simulated.

- **Number of dimensions**: Set the number of spatial dimensions (e.g., 2D or 3D).
  
- **Step length**: Set the step length for the Metropolis or Metropolis-Hastings algorithms, which controls the magnitude of each particle displacement.

- **Learning rate (for Gradient Descent programs)**: Specify the learning rate for gradient descent optimization of the variational parameters. A typical value is around \(1 \times 10^{-2}\).

- **Stopping criterion (for Gradient Descent programs)**: Define the stopping criterion based on energy convergence or the number of iterations. The process will stop once the change in energy is below a given threshold.

- **Oscillator frequency**: Set the frequency \(\omega\) of the harmonic oscillator potential that confines the particles.

- **Initial Wavefunction Variational parameters**: Modify the initial variational parameters (e.g., \(\alpha\), \(\beta_{ij}\)) for the wavefunction. These parameters will be optimized during the simulation to minimize the variational energy.

## Examples

### Example 1: Running a Serial Program for Non-Interacting Bosons with Gradient Descent

1. Copy the `main_NI_serial_GD.cpp` file from the `/VMC/MainFunctions` directory to the `/VMC` directory.
2. Compile and execute the program by running:

   ```bash
   ./execute


### Example 2: Running an MPI Program for Non-Interacting Bosons with Gradient Descent

1. Delete the previously copied `main_NI_serial_GD.cpp` file from the `/VMC` directory.
2. Copy the `9_boson_GD_MPI.cpp` file from the `/VMC/ResultGeneration/Boson` directory to the `/VMC` directory.
3. Compile and execute the MPI version of the program by running:

   ```bash
   ./execute_MPI


## License

This project is released into the public domain under the [Unlicense](https://unlicense.org).

Anyone is free to copy, modify, use, publish, compile, sell, or distribute this software, whether in source code or compiled binary form, for any purpose, commercial or non-commercial. There are no restrictions on how the software can be used.

In jurisdictions that recognize copyright, the authors dedicate any and all copyright interests in the software to the public domain, relinquishing all present and future rights under copyright law. This dedication is made for the benefit of the public and in perpetuity.

The software is provided "as is", without any warranty, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, or non-infringement. The authors are not liable for any claims, damages, or other liabilities arising from the use of the software.

For more details, please refer to [The Unlicense](https://unlicense.org).

## Acknowledgements

I would like to thank Øyvind Sigmundson Schøyen and Morten Ledum for their VMC solver template, which has been a valuable resource in developing this project.


## Contributing

We welcome contributions to improve this project. To contribute, follow these steps:

1. **Fork the repository** and create your own branch.
2. **Make your changes** or add new features.
3. **Write tests** to ensure your changes work as expected.
4. **Commit your changes** with clear and concise commit messages.
5. **Push to your fork** and submit a pull request to the `main` branch.

Please make sure your code follows the existing coding style and includes appropriate comments. Additionally, ensure that all tests pass before submitting a pull request.

For larger changes, it’s recommended to discuss the changes in an issue before starting your work.

We appreciate all contributions!
