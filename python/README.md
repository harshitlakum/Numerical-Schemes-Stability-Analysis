# Python Codes: Numerical Schemes and Stability Analysis

This folder contains Python scripts for:

1. **Fourth-order finite difference derivative (`fd_derivative_order.py`)**  
   Verifies the fourth-order accuracy of a central FD approximation to the first derivative. Produces error tables and convergence plots.

2. **Finite difference solver for 2D Poisson equation (`poisson2d_fd.py`)**  
   Solves the 2D Poisson equation with Dirichlet BCs using the standard 5-point stencil. Verifies convergence and plots numerical solution and error.

3. **Lax–Richtmyer stability demonstration (`lax_richtmyer_stability_demo.py`)**  
   Numerically demonstrates the Lax–Richtmyer stability bound for linear schemes.

## How to Run

Requires: numpy, matplotlib, scipy

Run any script with:

```bash
python script_name.py