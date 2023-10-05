#### Li Cycle Model ####

# import the deSolve library to numerically solve system
library(deSolve)

#### set model parameter values ####

# input and output fluxes
f_riv <- 1 # river
fht <- c(0.5, 1.3, 5, 10, 20) # hydrothermal vent
f_sr <- 0.6 # subduction reflux

# isotope ratios (isotopic composition) of input and output fluxes
d_riv <- 23 # riverine 
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
    
    # flux from clay sink defined here as it varies with f_ht 
    f_sc <- f_riv + f_ht + f_sr
    
    # rate of change
    ddelSW <- 1/N * (f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) + f_sc*(delSW - d_sc))
    
    # return the rate of change
    list(c(ddelSW))
  }) 
}

# set times parameter to output every 10 steps; for solver input
times <- seq(0, 5000000, 10000)

# f_ht = 0.5
f_ht <- fht[1]
parameters <- c(f_riv, f_ht, f_sr, d_riv, d_ht, d_sr, N)
out1 <- ode(y=state, times=times, func=Lithium, parms=parameters)
#out1 <- ode(y=state, times=times, func=Lithium, parms=parameters, method="rk4", hini=0.1)
out1 <- data.frame(out1)

# f_ht = 1.3
f_ht <- fht[2]
parameters <- c(f_riv, f_ht, f_sr, d_riv, d_ht, d_sr, N)
out2 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out2 <- data.frame(out2)

# f_ht = 5
f_ht <- fht[3]
parameters <- c(f_riv, f_ht, f_sr, d_riv, d_ht, d_sr, N)
out3 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out3 <- data.frame(out3)

# f_ht = 10
f_ht <- fht[4]
parameters <- c(f_riv, f_ht, f_sr, d_riv, d_ht, d_sr, N)
out4 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out4 <- data.frame(out4)

# f_ht = 20
f_ht <- fht[5]
parameters <- c(f_riv, f_ht, f_sr, d_riv, d_ht, d_sr, N)
out5 <- ode(y=state, times=times, func=Lithium, parms=parameters)
out5 <- data.frame(out5)

# plot resulting data with matplot function
plot(x=out1$time, y=out1$delSW, type="l", lwd=2, col="orange", ylim=c(24, 34))
lines(x=out1$time, y=out2$delSW, type="l", lwd=2, col="black")
lines(x=out1$time, y=out3$delSW, type="l", lwd=2, col="green")
lines(x=out1$time, y=out4$delSW, type="l", lwd=2, col="red")
lines(x=out1$time, y=out5$delSW, type="l", lwd=2, col="blue")
