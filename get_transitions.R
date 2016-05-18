# Function returning the probability of being in different states given previous time setPS
# and the neighborhood


######################################################
get_transitions = function(state0, EB, ET, EM, E1, E2, pars) 
{

    # COMPUTE THE LOGIT
    logit_alphab    = pars["ab0",1] + pars["ab1",1]*E1 + pars["ab2",1]*E2 + pars["ab3",1]*E1^2 + pars["ab4",1]*E2^2 + pars["ab5",1]*E1^3 + pars["ab6",1]*E2^3
    logit_alphat    = pars["at0",1] + pars["at1",1]*E1 + pars["at2",1]*E2 + pars["at3",1]*E1^2 + pars["at4",1]*E2^2 + pars["at5",1]*E1^3 + pars["at6",1]*E2^3
    logit_betab     = pars["bb0",1] + pars["bb1",1]*E1 + pars["bb2",1]*E2 + pars["bb3",1]*E1^2 + pars["bb4",1]*E2^2 + pars["bb5",1]*E1^3 + pars["bb6",1]*E2^3
    logit_betat     = pars["bt0",1] + pars["bt1",1]*E1 + pars["bt2",1]*E2 + pars["bt3",1]*E1^2 + pars["bt4",1]*E2^2 + pars["bt5",1]*E1^3 + pars["bt6",1]*E2^3
    logit_theta     = pars["th0",1] + pars["th1",1]*E1 + pars["th2",1]*E2 + pars["th3",1]*E1^2 + pars["th4",1]*E2^2 + pars["th5",1]*E1^3 + pars["th6",1]*E2^3
    logit_thetat    = pars["tt0",1] + pars["tt1",1]*E1 + pars["tt2",1]*E2 + pars["tt3",1]*E1^2 + pars["tt4",1]*E2^2 + pars["tt5",1]*E1^3 + pars["tt6",1]*E2^3
    logit_eps       = pars["e0",1]  + pars["e1",1]*E1 + pars["e2",1]*E2  + pars["e3",1]*E1^2 + pars["e4",1]*E2^2 + pars["e5",1]*E1^3 + pars["e6",1]*E2^3 

    #  Back transform into probabilities
    inv_logit = function(x) exp(x)/(1+exp(x))
    alphab = inv_logit(logit_alphab)
    alphat = inv_logit(logit_alphat)
    betab = inv_logit(logit_betab)
    betat = inv_logit(logit_betat)
    theta = inv_logit(logit_theta)
    thetat = inv_logit(logit_thetat)
    eps = inv_logit(logit_eps)
	phib = alphab*(EM + EB)*(1-alphat*(ET+EM))
	phit = alphat*(EM + ET)*(1-alphab*(EB+EM))
	phim = alphab*(EM + EB)*alphat*(EM + ET)
	
	p = matrix(0, nr = length(state0), nc = 4)

	# Compute transition probabilities		
	# From B
	p[state0=="B", 1] = (1 - eps - betat*(ET+EM)*(1-eps))[state0=="B"]
	p[state0=="B", 2] = 0
	p[state0=="B", 3] = (betat*(ET+EM)*(1-eps))[state0=="B"]
	p[state0=="B", 4] = eps[state0=="B"]

	# From T
	p[state0=="T", 1] = 0
	p[state0=="T", 2] = (1 - eps - betab*(EB+EM)*(1-eps))[state0=="T"]
	p[state0=="T", 3] = (betab*(EB+EM)*(1-eps))[state0=="T"]
	p[state0=="T", 4] = eps[state0=="T"]

	# From M
	p[state0=="M", 1] = (theta*(1-thetat)*(1-eps))[state0=="M"]
	p[state0=="M", 2] = (theta*thetat*(1-eps))[state0=="M"]
	p[state0=="M", 3] = ((1 - eps)*(1 - theta))[state0=="M"]
	p[state0=="M", 4] = eps[state0=="M"]

	# From R
	p[state0=="R", 1] = (alphab*(EM + EB)*(1-alphat*(ET+EM)))[state0=="R"]
	p[state0=="R", 2] = (alphat*(EM + ET)*(1-alphab*(EB+EM)))[state0=="R"]
	p[state0=="R", 3] = (alphab*(EM + EB)*alphat*(EM + ET))[state0=="R"]
	p[state0=="R", 4] = (1 - rowSums(p[state0=="R",1:3]))
	
	return(p)
}