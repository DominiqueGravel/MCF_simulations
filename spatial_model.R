# pres: vector of states for each class
# E1: annual average temperature (scaled)
# E2: annual preciptation (scaled)
# XY: X and Y coordinates for each cell
# pars: parameters for the transitions
# nsteps: number of steps (5yr) to run the model

spatial_model = function(pres, E1, E2, XY, pars, nsteps) {

	# Prepare the object to store the time series
	series = matrix(nr = nsteps, nc = 4)

	# Compute the connectivity matrix
	N = length(pres)
	distMat = as.matrix(dist(XY, method = "euclidean", upper = T, diag = T))
	ConMat = matrix(0, nr=N,nc=N)
	ConMat[distMat<1.5] = 1
	diag(ConMat) = 0

	# MAIN LOOP OVER TIME
	for(step in 1:nsteps) {

		# Transform states into binary presence-absence
		presB = numeric(N)
		presB[pres == "B"] = 1
		presT = numeric(N)
		presT[pres == "T"] = 1
		presM = numeric(N)
		presM[pres == "M"] = 1		
	
		# Compute occupancy in the neighbourhood
		ConB = (ConMat%*%presB)/8
		ConT = (ConMat%*%presT)/8
		ConM = (ConMat%*%presM)/8

		# Compute the probability of being in each state at t1
		transitions = get_transitions(pres, ConB, ConT, ConM, E1, E2, pars) 
		
		# Draw the identify of the state at t1
		pres = apply(transitions, 1, draw)
	
		# Save regional abundnace of each state
		series[step,] = c(sum(pres%in%"B"),sum(pres%in%"T"),sum(pres%in%"M"),sum(pres%in%"R"))/N

		cat(step,'\n')
		}

	return(list(map = pres, series=series))
}	

	
################################	
draw = function(p) {
	states = c("B","T","M","R")
	draw = rmultinom(n=1,size=1,prob=p)
	states[which(draw==1)]
}


