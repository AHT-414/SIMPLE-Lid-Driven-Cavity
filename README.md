# SIMPLE-Lid-Driven-Cavity
2D cavity flow using the SIMPLE Algorithm

This 2D Lid Driven Cavity flow is written in FORTRAN and uses the SIMPLE algorith to numerically solve the NS eqautions.

The main file is cavityFlow.f95, other files included are modules and subroutines. The solvePressureMatrix subroutine includes functions
from LAPACK (http://www.netlib.org/lapack/) to solve the system of equations. The program cannot be compiled without installing the package.

A makefile is also included to create the executable and use the correct flags to include the LAPACK functions.
