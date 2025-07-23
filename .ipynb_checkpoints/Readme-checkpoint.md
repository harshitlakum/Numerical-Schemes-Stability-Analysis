# Numerical Schemes & Stability Analysis

This repository provides a collection of Python and MATLAB scripts for the study and demonstration of core numerical analysis concepts, including finite difference schemes, convergence analysis, and stability theorems.

## Contents

**Python and MATLAB codes are provided for:**
- **High-order finite difference derivative approximations** (and numerical order verification)
- **Finite difference method for the 2D Poisson equation** (with Dirichlet boundary conditions and error/convergence analysis)
- **Lax–Richtmyer stability theorem** (numerical demonstration for linear time-stepping schemes)

---

## Directory Structure

Numerical-Schemes-Stability-Analysis/
│
├── python/
│   ├── fd_derivative_order.py
│   ├── poisson2d_fd.py
│   ├── lax_richtmyer_stability_demo.py
│   ├── README.md
│
├── matlab/
│   ├── fd_derivative_order.m
│   ├── poisson2d_fd.m
│   ├── lax_richtmyer_stability_demo.m
│   ├── README.md
│
├── LICENSE
└── README.md        ← (this file)

---

## Scripts Overview

| Script                         | Description                                 | Language |
|--------------------------------|---------------------------------------------|----------|
| fd_derivative_order            | 4th-order finite difference derivative and accuracy/convergence verification | Python / MATLAB |
| poisson2d_fd                   | 2D Poisson equation (finite difference, Dirichlet BC, convergence) | Python / MATLAB |
| lax_richtmyer_stability_demo   | Demonstrates the Lax–Richtmyer stability bound numerically | Python / MATLAB |

---

## How to Use

### Prerequisites

- **Python scripts:**  
  Require `numpy`, `matplotlib`, and `scipy`. Install with  
  ```bash
  pip install numpy matplotlib scipy

	•	MATLAB scripts:
Require basic MATLAB (no toolboxes needed beyond core).

⸻

Running the Scripts

Python

From the python/ directory, run for example:

python fd_derivative_order.py
python poisson2d_fd.py
python lax_richtmyer_stability_demo.py

MATLAB

From the matlab/ directory or with that as your MATLAB working folder:

fd_derivative_order
poisson2d_fd
lax_richtmyer_stability_demo


⸻

What You’ll See
	•	fd_derivative_order: Prints error and observed convergence order for the 4th-order finite difference formula; produces a log-log error plot.
	•	poisson2d_fd: Solves the 2D Poisson problem on a unit square, plots the solution and error, and prints/plots the error vs. mesh size (showing second-order convergence).
	•	lax_richtmyer_stability_demo: Evolves a linear system and verifies that the Lax–Richtmyer stability bound holds numerically.
