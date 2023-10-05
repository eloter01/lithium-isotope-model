#### Li Cycle Model ####

# import the deSolve library to numerically solve system
library(deSolve)

#### set model parameter values ####

# input and output fluxes
f_riv <- 15 # river
f_ht <- 5 # hydrothermal vent
f_sr <- 0.6 # subduction reflux
f_sc <- 20.6 # Li sink into clays; clay formation during seafloor alteration and sediment diagenesis 

# isotope ratios (isotopic composition) of input and output fluxes
d_riv <- 2 # riverine 
d_ht <- 8 # hydrothermal vent
d_sr <- 15 # subduction reflux 

# constant seawater total Li
N <- 3.6e6

# set parameters to vector for solver input
parameters <- c(f_riv, f_ht, f_sr, f_sc, d_riv, d_ht, d_sr, N)

# state variables for solver input
state <- c(delSW = 31)

# model equation for solver input
Lithium <- function(t, state, parameters) {
  with(as.list(c(t, state, parameters)),{
    
    d_sc <- delSW - 16
    
    # rate of change
    ddelSW <- 1/N * (f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) + f_sc*(delSW - d_sc))
    
    # return the rate of change
    list(c(ddelSW))
    }) 
}

# set times parameter to output every 10 steps; for solver input
times <- seq(0, 5000000, 10000)

# run solver and save output
out <- ode(y=state, times=times, func=Lithium, parms=parameters)

tail(out)

# plot resulting data with matplot function
plot(x=out[,1], y=out[,2], type="l", lwd=2, col="blue", xlim=c(0, 5000000), ylim=c(15, 35))
