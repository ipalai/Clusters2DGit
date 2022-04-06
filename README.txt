This repository contains code to repeat the simulations described in the following paper:

 I. Palaia and A. Saric
 Controlling cluster size in 2D phase-separating binary mixtures with specific interactions
 https://doi.org/10.1101/2022.02.10.479877


--------------
 DESCRIPTION
--------------
The code creates LAMMPS input scripts. The system contains 2 types of patchy particles (A and B). They have finite valence q_A = q_B = q and can form specific heterotypic bonds (A with B). 

In StaticClusters, no isotropic attraction between particles is present, so that clusters are solid.
In FluidCLusters, an isotropic attraction (A-B, A-A and B-B) is added, so that clusters are fluid.

See the Methods section of the paper for more details.


--------------
 INSTRUCTIONS
--------------

Run the code like this:
 python3 make_simulation.py
This will create two files. The Input*.in file is the LAMMPS input file, while the Config*.dat file contains an initial configuration for the system with isolated particles placed at random on a grid.

The previous instruction will use default parameters. Parameters can be changed directly from the command line like this:
 python3 make_simulation.py -nA 200 -nB 50 -eps_ll 10 -density_A_target 0.03 -real 100 -qA 5 -qB 5 -runsteps 50000 -dumpstep 500 -ea 2 -ra 2

To know all available options run:
 python3 make_simulation.py --help


Note that parameters are named in a slightly different manner in the code and in the paper. Particular care is needed when looking at lengths, as the diameter of a particle is taken equal to 2 in simulations while it is taken equal to 1 sigma in the paper. Here is an equivalence table (the last 2 lines are only relevant for fluid clusters).
				SIMULATIONS		PAPER
particle diameter		sp = 2.0		sigma
inner patch radius		sl = 0.1		r_in = 0.05 sigma
outer patch radius		sl+rl = 0.3		r_out = 0.15 sigma
range of patch attraction	rl = 0.2		r_out-r_in = 0.1 sigma
mass				1			m
density				0.0075			0.03 sigma^(-2)
outer isotropic attr. radius	sp+ra (= 4 e.g.)	R_a (= 2.0 sigma e.g.)
range of isotropic attraction	ra (= 2 e.g.)		R_a - sigma ( = 1.0 sigma e.g.)

