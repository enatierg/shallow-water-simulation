# One-Dimensional Shallow Water Model

This repository contains the implementation of a **one-dimensional Shallow Water model (SWM)** using **FORTRAN95** with double precision arithmetic. Notably, part of the task was to improve certain parts of the code, so I can't claim ownership of the full code. The parts done by me are indicated by the comments.

The project includes:  

* **Matrix setup and linear solve**: Construction of the finite difference matrix `Dx` and Helmholtz matrix `H`, along with the `helmholtz_solve` subroutine to solve `H * eta = r`. Includes testing of the solver with random input vectors.  
* **Time-stepping algorithm**: Implementation of the theta-method for the shallow water equations, using subroutine `timestepping` to update velocity and height perturbations over a given number of timesteps.  
* **Initial conditions**: Subroutine `initial_condition` sets up the initial velocity and height profiles.  
* **Output and visualization**: Subroutine `save_fields` writes computed fields to text files, with an accompanying `visualise.py` script to plot results.  
* **Field comparison**: Function `field_delta` calculates the difference between two sets of fields to quantify numerical accuracy.  

The code is structured for **clarity and efficiency**, exploiting properties of banded and symmetric matrices and using BLAS/LAPACK routines where appropriate. A Makefile is included for compilation.   
