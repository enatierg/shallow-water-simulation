# Two-Dimensional Shallow Water Model

This repository contains a parallel implementation of a **two-dimensional Shallow Water model** developed as part of a scientific computing course. The project demonstrates advanced numerical methods for solving partial differential equations using **FORTRAN95** with MPI for distributed computation.

## Technical Overview

This implementation was developed as part of a high-performance computing curriculum focusing on:
- Parallel algorithms for scientific computing
- Distributed memory programming with MPI
- Numerical methods for fluid dynamics
- Conservation properties of numerical schemes
- Performance optimisation techniques

## Key Variables

- **m**: Number of grid points in one direction
- **kmax**: Number of time steps
- **nu**: Courant number ν = cg∆t/∆x
- **Dx**: Sparse real distributed n × n Matrix for the difference operator Dx
- **Dy**: Sparse real distributed n × n Matrix for the difference operator Dy
- **u**: Real distributed Vector of length n of the discretised horizontal velocity u
- **v**: Real distributed Vector of length n of the discretised vertical velocity v
- **eta**: Real distributed Vector of length n of the height perturbation η
- **rank**: MPI rank of the current process
- **nproc**: Total number of processes

## Project Structure

```
├── sw_main.f90          # Main program
├── header.f90           # Type definitions and module declarations
├── matmult.f90          # Matrix-vector multiplication operations
├── eta_peak.f90         # Initial condition setup
├── average.f90          # Parallel averaging across processes
├── leapfrog.f90         # Leapfrog time integration algorithm
├── create_matrices.f90  # Finite difference matrix setup
├── sparsegather.f90     # MPI communication routines
├── save_field.f90       # Data output functionality
├── test_leapfrog.f90    # Validation and testing suite
├── visualise.py         # Python visualisation script
└── Makefile             # Build configuration
```

## Key Components

* **Matrix Setup**  
  Implementation of finite difference matrices `Dx` and `Dy` for spatial derivatives using periodic boundary conditions in both dimensions.

* **Parallel Communication**  
  MPI-based communication routines (`sparsegather`) for exchanging boundary data between processes using an alternating pattern to avoid deadlocks.

* **Time Integration**  
  Implementation of the leapfrog algorithm with alternating κ values for stable time integration of velocity and height fields.

* **Conservation Checking**  
  Verification of mass conservation through the `average` subroutine that computes global averages across all MPI processes.

* **Testing Framework**  
  Comprehensive test suite (`test_leapfrog`) to validate both the time integration algorithm and averaging functionality.

## Algorithm Overview

The implementation follows a leapfrog algorithm with:
- Initial conditions featuring a cosine-shaped peak in height field (η)
- Zero initial velocity fields (u, v)
- Conservation checking of the discrete integral I(k) = (Δx)²Σηᵢⱼ(k)
- Periodic boundary conditions in both spatial dimensions

## Usage

Compile the code using the provided Makefile:
```bash
make
```

Main executables:
- `sw_main`: Primary simulation program
- `test_leapfrog`: Validation and testing programme

## Technical Features

- **Parallel Design**: Domain decomposition with MPI processes
- **Matrix Storage**: Compressed row storage format for efficient sparse matrix operations
- **Boundary Handling**: Special treatment for periodic boundaries and inter-process communication
- **Performance Timing**: MPI-based timing of computational performance

The implementation maintains the conservation properties of the continuous equations in the discrete formulation and efficiently handles the sparse matrix-vector operations required for the finite difference discretisation.
