Simulation code for the MCF summer school
============

## Getting started

#### 1. Download repository

	git clone git@github.com:DominiqueGravel/MCF_simulations.git
	cd MCF_simulations

#### 2. Run the simulation model

You have to open the "run_simulations.r" script and play with the code in it. There is already an example.

#### 3. How to tweak the model

The rules determining the transitions are provide in the file "get_transistions.r". The parameters are the result of a fairly complex statistical analysis and it is not recommended to play with them. It is however feasible to manually modify the probabilities to force scenarios of forest management, such as plantation (changing the R->T and R->B transitions), selective logging (changing the transitions M->T and M->B), along with clearcuts (changing the M->R, T->R and B->R transitions).

It is also possible to implement assisted migration by modifying the connectivity matrix. This could be done for instance by changing the threshold distance determining if two cells of the lattice are connected by dispersal. It could also be possible to implement density-independent dispersal, but this would require playing with the probabilities (high-level tweaking).
