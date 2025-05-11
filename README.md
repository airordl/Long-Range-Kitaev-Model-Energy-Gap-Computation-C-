Long-Range Kitaev Model — Energy Gap Computation (C++)

This C++ project implements a solver for computing the energy spectrum and energy gap of the long-range Kitaev model on a 2D lattice.
It uses exact diagonalization of the Hamiltonian, which includes long-range interactions extending beyond nearest neighbors. 
The solver is modular, efficient, and written with custom matrix manipulation routines for portability and performance.

Overview:
- Solves the Kitaev model with long-range spin-spin couplings
- Constructs and diagonalizes the full Hamiltonian matrix
- Computes the many-body energy gap for finite lattice configurations
- Implements core matrix algebra internally (no external libraries)

Files:
- `kitaevgap.cc` — Main program to set up and solve the long-range Kitaev Hamiltonian
- `kitaev.h` — Header file defining model parameters and Hamiltonian construction
- `matx.cc` — Implementation of basic matrix operations
- `matx.h` — Matrix class and helper utilities (e.g., dot products, diagonalization)

How to Build:
Compile using a standard C++ compiler

