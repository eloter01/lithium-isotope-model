#### Li Cycle Model ####

# import the deSolve library to numerically solve system
library(deSolve)


#### set model parameter values ####

# input and output fluxes
f_riv <- 1 # river
f_ht <- 13 # hydrothermal vent
f_sr <- 0.6 # subduction reflux

# isotope ratios (isotopic composition) of input and output fluxes
d_riv <- 23 # riverine 
d_ht <- 8 # hydrothermal vent
d_sr <- 15 # subduction reflux 



## uncomment this section for deSolve model

# input <- f_riv + f_ht + f_sr
# input2 <- f_riv + f_ht + f_sr
# f_sc <- 2.9
# output2 <- f_sc
# N0 <- 3.4e6
# 
# # state variables for solver input
# state <- c(delSW = 31, N = 3.4e6)
# 
# # model equation for solver input
# Lithium <- function(t, state, parameters) {
#   with(as.list(c(t, state, parameters)),{
#     
#       # isotope composition of clay sink defined here as it varies with delSW
#       d_sc <- delSW - 16
# 
#       f_sc <- 2.9 * ((input2 - output2) + 3.4e6)/3.4e6
#       input2 <- input2 + input
#       output <- f_sc
#       output2 <- output2 + f_sc
# 
#       # rate of change
#       ddelSW <- 1/N * (f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) + f_sc*(delSW - d_sc))
#       dN <- input - output
#       
#     # return the rate of change
#     list(c(ddelSW, dN))
#   }) 
# }

# set times parameter to output every 10000 steps; for solver input
times <- seq(0, 5000000, 10000)



###############
## comment this section when using deSolve model

input <- f_riv + f_ht + f_sr
input2 <- f_riv + f_ht + f_sr
f_sc <- 2.9
output2 <- f_sc
delSW <- 31
dN <- 3.4e6
N0 <- 3.4e6

for (i in 1:400){
  # isotope composition of clay sink defined here as it varies with delSW
  d_sc <- delSW - 16

  f_sc <- 2.9 * ((input2 - output2)*10000 + 3.4e6)/3.4e6
  input2 <- input2 + input
  output <- f_sc
  output2 <- output2 + f_sc

  # rate of change
  delSW <- delSW + 10000/N0 * (f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) + f_sc*(delSW - d_sc))
  dN <- dN + input*10000 - output*10000
}
################



## uncomment this section for deSolve model

# parameters <- c(f_riv, f_ht, f_sr, d_riv, d_ht, d_sr, N0)
# out1 <- ode(y=state, times=times, func=Lithium, parms=parameters, method="rk4")
# out1 <- data.frame(out1)
# 
# plot(x=out1$time, y=out1$delSW, type="l", lwd=2, col="orange")
