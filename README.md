# AtomicKohnSham.jl

[![Build Status](https://github.com/Theozeud/AtomicKohnSham/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Theozeud/AtomicKohnSham/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Theozeud/AtomicKohnSham/branch/main/graph/badge.svg)](https://codecov.io/gh/Theozeud/AtomicKohnSham)
![License](https://img.shields.io/badge/license-MIT-blue.svg)

**AtomicKohnSham.jl** is a Julia package designed to **compute the ground state of isolated atoms and ions** within the framework of **Extended Kohn-Sham models**. Originally developed to investigate the existence of negative ions, the code also serves as a versatile tool for testing new density functionals and exploring numerical precision issues.

These models require solving a nonlinear eigenvalue partial differential equation (PDE) in three dimensions. The nonlinearity arises from the Kohn–Sham potential, which depends on a finite number of eigenfunctions. As is standard in computational chemistry, a fixed-point iterative scheme (Self-Consistent Field or SCF) is used, alternating with the solution of the eigenvalue problem.

Assuming spherical symmetry, the eigenvalue problem can be reduced to a family of radial 1D equations. This package implements a **high-precision finite element method (FEM)** for solving these equations, using **polynomials of degree up to 20** on exponential meshes. This approach allows for highly accurate results with relatively few mesh points, often outperforming low-order methods on refined meshes.

The code supports both **double (Float64)** and **quadruple (Double64)** floating-point precision. This level of precision is crucial for resolving subtle numerical questions that remain open in the literature. For example, in certain Extended Kohn–Sham models, it is still unclear whether the energy of the outermost orbital in some atoms is exactly zero or simply very close to zero but negative—a distinction that this package is designed to investigate numerically.

## Features
Currently supported functionalities:

- **SCF Algorithms**:
  - Constant Damping Algorithm (CDA)
  - Optimal Damping Algorithm (ODA)
  - Quadratic damping schemes
    
- **Symmetry**:
  - Spherical symmetry reduction to radial equations
    
- **Finite Element Basis**:
  - Integrated Legendre polynomials (up to order 20)
    
- **Exchange–Correlation Functionals**:
    - LDA, LSDA (via Libxc.jl)
      
- **Mesh types**:
    - Linear, geometric, exponential
      
- **Floating-point precision**:
    - Float64, Double64 (quadruple precision via DoubleFloats.jl)

## Recommendations
Based on our experiments, the following setup provides excellent results:
- Use an exponential mesh with a parameter between 1 and 2.
- Combine it with the Integrated Legendre polynomial basis of order up to 20.
  
This combination offers a good balance between numerical accuracy and computational efficiency.

## Notes :
- Currently, only one FEM basis has been implemented. However, adding a new one is straightforward, as all FEM matrix assemblies are handled automatically


- Higer precision should works but this has not been thoroughly tested.

## Example

Here is one example on the hydrogen atom with the Reduced-Hartree Fock model. 
```julia
problem = AtomProblem(;
                T               = Float64,                    # Data type for computations
                model           = RHF(;z=1, N=1),             # Model : Reduced-Hartree Fock
                                                              # with nuclear chare z=1  and N=1 electrons
                                                              # (= Hydrogen)
                name            = "Hydrogen",                 # A custom name that you chan choose
                
                alg             = ODA(0.4),                   # SCF algorithm : Optimal Dampling with initial
                                                              # parameter equal to 0.4
                scftol          = 1e-11,                      # Tolerance for the scf procedure
                maxiter         = 60,                         # Maximum number of SCF iterations :
                                                              # If this number is reached, the algorithm stops
                                                              # regardless of convergence 
                degen_tol       = 1e-2,                       # Tolerance to detect degeneracy between orbital energies

                lh              = 0,                          # Angular momentum cutoff (ℓ ≤ lh)
                                                              # only the s-orbital is needed for the hydrogen atom
                Rmax            = 100,                        # Radial domain cutoff
                Nmesh           = 10,                         # Number of points of the mesh
                typemesh        = expmesh,                    # The mesh used is an exponential mesh
                optsmesh        = (s = 1.5,),                 # Mesh parameters: here, s = 1.5
                typebasis       = P1IntLegendreGenerator,     # The FEM Basis is composed of the integrated legendre 
                                                              # polynomials with the P1 elements
                optsbasis       = (ordermax = 20,),           # Polynomials up to order 20 are used

                verbose         = 0)                          # Verbosity level: 0 = silent, 3 = maximum verbosity
```
Additional options are available but have been omitted here for simplicity.
