# UniProjects
University Assignments covering numerical methods for solving differential equations, and parallel computing method in both OpenMP and OpenMPI.

"Scientific Computing" labelled codes cover my assignment solutions for my third year scientific computing module focusing on numerical methods for solving ordinary and partial differential equations subject to boundary conditions. 

Of particular interest is "Tokamak_HeatFLow.c" which models the heat flow around the torodial shell of a tokamak in both a steady-state (i.e. at the late time, after the temperature distribution has relaxed) and in the time-evolution mode for the initial condition of T=0.The program reads in a series of input parameters from 'input.txt' and reads in a series of position-variant functions from 'coefficients.txt' and outputs the temperature profile for all grid points in one single column to the file 'output.txt'. For the case of the steady state, the PDE is  solved through using a banded system of linear equations, discretised using finite difference methods, and then solved through using LAPACKE's dgbsv.f. While for the case of the time-evolution mode, an implicit solver is employed in combination with a finite difference method. These routines are adapted from the lectures notes of Prof. McMillan's PX390 module. 


"HPC:*" labelled codes contain my assigment solutions for my final year high performance and parallel computing modules, which covered a range of methods including OpenMP, OpenMPI and CUDA. Codes labelled "HPC" contain some of the neccesary files to generate random numbers and constructs pngs of solutions.

"ofcmodel.c" contains a hybrid OpenMP-MPI parallelisation of the Olami-Feder-Christensen Model: used to model earthquakes through self-organised criticality. The OFC model attempts to model the stresses along a fault line, with distance along the fault and depth defining a 2D grid of sites. The stresses on these sites evolve according to simple rules with just one tunable parameter, α plus the grid sizes. These rules are enough to generate SOC across a range of values of α.
