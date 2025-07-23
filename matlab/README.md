# MATLAB Codes: Numerical Schemes and Stability Analysis

This folder contains MATLAB scripts for:

1. **Fourth-order finite difference derivative (`fd_derivative_order.m`)**  
   Verifies the fourth-order accuracy of a central FD formula for the first derivative. Prints error table and plots convergence.

2. **Finite difference solver for 2D Poisson equation (`poisson2d_fd.m`)**  
   Solves the 2D Poisson equation with Dirichlet BCs (5-point stencil). Shows numerical solution, error, and convergence study.

3. **Lax–Richtmyer stability demonstration (`lax_richtmyer_stability_demo.m`)**  
   Numerically demonstrates the Lax–Richtmyer stability theorem for a linear evolution scheme.

## How to Run

Run any script by name in MATLAB (from this directory), e.g.:

```matlab
fd_derivative_order
poisson2d_fd
lax_richtmyer_stability_demo