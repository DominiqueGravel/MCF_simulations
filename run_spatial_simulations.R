rm(list=ls())

# Parameters
pars = read.table("pars.txt",row.names=1)

# Functions
source("get_transitions.R")
source("spatial_model.R")

# SETTING THE LATTICE
Xlim = 250
Ylim = 50
N = Xlim*Ylim
X = c(1:Xlim)
Y = c(1:Ylim)
XY = expand.grid(X,Y)

# SCALED ENVIRONMENTAL VARIABLES
load("scale_info.Robj")
tmin = (-3 - vars.means[1])/vars.sd[1]
tmax = (10 - vars.means[1])/vars.sd[1]
E1 = rep(seq(tmax,tmin, length=Xlim), times = Ylim)
E2 = numeric(N)

# CONVERT INTO REAL VALUES
E1_real = E1*vars.sd[1] + vars.means[1]
E2_real = E2*vars.sd[2] + vars.means[2]

# INITIAL CONDITIONS
initial_probs = matrix(0.25, nr = N, nc = 4)
pres = apply(initial_probs, 1, draw)

# RUN THE MODEL
large_run = spatial_model(pres = pres, E1 = E1, E2 = E2, XY=XY, pars=pars, nsteps = 100)

# PLOT THE RESULTS
dev.new(width = 5, height = 15)
X = sort(unique(XY[,1]))
Y = sort(unique(XY[,2]))
pres_int = as.factor(large_run$map)
n = 1
Z = matrix(nr = length(X), nc = length(Y))
for(x in 1:length(X)) 
	for(y in 1:length(Y)) 	
		Z[x,y] = pres_int[XY[,1] == X[x] & XY[,2]==Y[y]]	
	

layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(2,12))
par(mar=c(0,0,2,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title(2015,cex=2)
legend("center",legend = c("Boreal","Mixte","Regen","Temperate"),fill = c("darkcyan","palegreen3","white","orange"),bty = "n",cex = 1.5)
par(mar=c(5,5,0,5))
image(Y,X,t(Z),cex.axis = 1.5, cex.lab = 1.25, col = c("darkcyan","palegreen3","white","orange"))


