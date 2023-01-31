# FLDCRV-GRID
A program utilizing the core functionality of the FLDCRV program written by J.D. Landstreet, 1981, which computes synthetic longitudinal field curves in a grid and compares them against input (real) data.

# FLDCRV
FLDCRV is a program written by J.D. Landstreet in 1981. It originally computed synthetic longitudinal field curves over a fixed, equally spaced grid from -0.5 to 0.5 in rotational phase.

# My adjustments
I've re-written the program, utilizing only the core functionality provided by the unfld and magfld subroutines (which have not been provided in this repo as it's not my code to share).
My program imports a set of real, phased longitudinal field observations and computes synthetic longitudinal field profiles according to a grid in inclination, beta angle, and dipolar field strength, determining the best-fit parameters through chi-squared minimization.
