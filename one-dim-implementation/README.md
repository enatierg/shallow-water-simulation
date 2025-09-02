# One-Dimensional Shallow Water Model

This folder contains the implementation of a **one-dimensional Shallow Water model (SWM)** using **FORTRAN95** with double precision arithmetic. Not all parts of the code were implemented by me; the parts completed by me are indicated below.  

## Project Components

* **Matrix setup and linear solve**  
  Construction of the finite difference matrix `Dx` and Helmholtz matrix `H`.  
  **Implemented by me**: the `helmholtz_solve` subroutine to solve `H * eta = r` and tests to verify its correctness with random input vectors.  

* **Time-stepping algorithm**  
  **Implemented by me**: the main theta-method loop in the subroutine `timestepping`, which updates velocity (`u`) and height perturbation (`eta`) over a given number of timesteps.  

* **Initial conditions**  
  Subroutine `initial_condition` sets up the initial velocity and height profiles. Adapted by me for integration with the timestepping routine.  

* **Output and visualization**  
  Subroutine `save_fields` writes computed fields to text files. Accompanied by `visualise.py` to plot results.  

* **Field comparison**  
  Function `field_delta` calculates the difference between two sets of fields to quantify numerical accuracy.  

## Notes

The code is structured for **clarity and efficiency**, exploiting properties of banded and symmetric matrices and using BLAS/LAPACK routines where appropriate. A Makefile is included for compilation.
