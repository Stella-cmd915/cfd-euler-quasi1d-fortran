# Quasi-1D Euler Solver (Fortran)

A finite-volume solver for compressible flow in a converging–diverging duct.

## Tech stack
- Fortran 90

## Features
- Roe and Steger–Warming flux-vector splitting
- RK4 time integration
- 1st/2nd-order MUSCL with Van Leer limiter
- Input cases in `inputs/`

## Getting started
```bash
gfortran -O2 src/euler_q1d.f90 -o euler_q1d
./euler_q1d < inputs/inp_roe1.txt
```

---
Part of the **SIGMA LAB** project portfolio.
