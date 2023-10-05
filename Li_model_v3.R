# Solving Li box model by 4th order Runge-Kutta method


#########      Main Equation        #########

# One box model

#dNdt <- f_riv + f_ht + f_sr - f_sc
#ddeltaSWdt <- 1/N*(f_riv*(d_riv - deltaSW) + f_ht*(d_ht - deltaSW) + f_sr*(d_sr - deltaSW) + f_sc*(d_sc - deltaSW))


#########      Define constants        #########

# input and output fluxes
f_riv = 1 # river
f_ht = 1.3 # hydrothermal vent
f_sr = 0.6 # subduction reflux
f_sc = f_riv + f_ht + f_sr # Li sink into clays; clay formation during seafloor alteration and sediment diagenesis 

# isotope ratios (isotopic composition) of input and output fluxes
#d_riv = 23 # riverine 
d_ht = 8 # hydrothermal vent
d_sr = 15 # subduction reflux 
d_sc = 15 # Li sink into clays


#########      Define two main functions       #########
# t is timestep, M is sewater Li reservoir, and deltaSW is seawater Li isotopic balance
fM <- function(t, M, deltaSW) 
{
  f_riv + f_ht + f_sr - f_sc
}

fdeltaSW <- function(t, M, deltaSW, driv) 
{
  1/M*(f_riv*(driv - deltaSW) + f_ht*(d_ht - deltaSW) + f_sr*(d_sr - deltaSW) - f_sc*(d_sc - deltaSW))
}

# d_riv is a function of time
d_riv_func <- function(t){
  if (t<1000000){
    d_riv <- 24.1
  } else if (t>=1000000 && t<2000000){
    d_riv <- 15
  } else if (t>=2000000 && t<3000000){
    d_riv <- 8
  } else if (t>=3000000 && t<4000000){
    d_riv <- 5
  } else if (t>=4000000){
    d_riv <- 2
  }
  
  return(d_riv)
}

#define step size
dt=10000
t_total=5000000
t_start = 0
N=(t_total/dt)
t = seq(1:N)

#initial conditions
M = seq(1:N) #initial total Ca in seawater
deltaSW = seq(1:N) #initial d44Ca of seawater
t[1]=t_start
M[1]=3.4e6
deltaSW[1]=31


###########    4th order Runge-Kutta Solver ###########
rivs <- c(0)
for (i in 1:N)
{
  #update time
  t[i+1] <- t[i] + dt
  
  # calculate d_riv as function of time
  driv <- d_riv_func(t[i])
  
  
  # vector to see that driv is changing correctly according to rules
  rivs <- c(rivs, driv)
  
  #update M & deltaSW
  k1M <- fM(t[i], M[i], deltaSW[i])
  k1deltaSW <- fdeltaSW(t[i], M[i], deltaSW[i], driv)
  
  k2M <- fM(t[i]+dt/2, M[i]+dt/2*k1M, deltaSW[i]+dt/2*k1deltaSW)
  k2deltaSW <- fdeltaSW(t[i]+dt/2, M[i]+dt/2*k1M, deltaSW[i]+dt/2*k1deltaSW, driv)
  
  k3M <- fM(t[i]+dt/2, M[i]+dt/2*k2M, deltaSW[i]+dt/2*k2deltaSW)
  k3deltaSW <- fdeltaSW(t[i]+dt/2, M[i]+dt/2*k2M, deltaSW[i]+dt/2*k2deltaSW, driv)
  
  k4M <- fM(t[i]+dt, M[i]+dt*k3M, deltaSW[i]+dt*k3deltaSW)
  k4deltaSW <- fdeltaSW(t[i]+dt, M[i]+dt*k3M, deltaSW[i]+dt*k3deltaSW, driv)
  
  M[i+1] <- M[i]+dt/6*(k1M + 2*k2M + 2*k3M + k4M)
  deltaSW[i+1] <- deltaSW[i]+dt/6*(k1deltaSW + 2*k2deltaSW + 2*k3deltaSW + k4deltaSW)
}


##########  plot numerical results ############

#plot "Mass of Li in the seawater"
plot(t, M, type="l")

#plot " Li isotopic balance of seawater"
plot(t, deltaSW, type="l")


