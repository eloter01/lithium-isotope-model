#### Li Cycle Model ####

# import the deSolve library to numerically solve system
library(deSolve)

#### set model parameter values ####

# input and output fluxes
f_riv <- 1 # river
f_ht <- 1.3 # hydrothermal vent
f_sr <- 0.6 # subduction reflux
f_sc <- f_riv + f_ht + f_sr # Li sink into clays; clay formation during seafloor alteration and sediment diagenesis 

# isotope ratios (isotopic composition) of input and output fluxes
driv <- c(20, 10, 5, 2) # riverine 
d_ht <- 8 # hydrothermal vent
d_sr <- 15 # subduction reflux 

# constant seawater total Li
N <- 3.6e6

# state variables for solver input
state <- c(delSW = 31)

# model equation for solver input
Lithium <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # isotope composition of clay sink defined here as it varies with delSW
    d_sc <- delSW - 16
    
    # rate of change
    ddelSW <- 1/N * (f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) + f_sc*(delSW - d_sc))
    
    # return the rate of change
    list(c(ddelSW))
  }) 
}

# set times parameter to output every 10 steps; for solver input
times <- seq(0, 5000000, 10000)

# d_riv = 20
d_riv <- driv[1]
parameters <- c(f_riv, f_ht, f_sr, f_sc, d_riv, d_ht, d_sr, N)
out1 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out1 <- data.frame(out1)

# d_riv = 10
d_riv <- driv[2]
parameters <- c(f_riv, f_ht, f_sr, f_sc, d_riv, d_ht, d_sr, N)
out2 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out2 <- data.frame(out2)

# d_riv = 5
d_riv <- driv[3]
parameters <- c(f_riv, f_ht, f_sr, f_sc, d_riv, d_ht, d_sr, N)
out3 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out3 <- data.frame(out3)

# d_riv = 2
d_riv <- driv[4]
parameters <- c(f_riv, f_ht, f_sr, f_sc, d_riv, d_ht, d_sr, N)
out4 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out4 <- data.frame(out4)

# plot resulting data with matplot function
plot(x=out1$time, y=out1$delSW, type="l", lwd=2, col="orange", ylim=c(22, 32))
lines(x=out1$time, y=out2$delSW, type="l", lwd=2, col="green")
lines(x=out1$time, y=out3$delSW, type="l", lwd=2, col="red")
lines(x=out1$time, y=out4$delSW, type="l", lwd=2, col="blue")
