# Solving Li box model by 4th order Runge-Kutta method


#########      Main Equation        #########

# One box model

#dNdt <- f_riv + f_ht + f_sr - f_sc
#ddelSWdt <- 1/N*(f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) + f_sc*(d_sc - delSW))



#########      Define constants and parameters       #########


# define time domain for 4th Order Runge-Kutta solver
dt <- 10000 # time step
t_total <- 5000000 # total num yrs
t_start <- 0 # start time
t <- seq(t_start, t_total, by=dt) # time sequence
N_step <- length(t) - 1 # num of time steps for solver

# input and output fluxes
friv = c(2.5, 5, 10, 50) # river
f_ht = 1.3 # hydrothermal vent
f_sr = 0.6 # subduction reflux

# isotope ratios (isotopic composition) of input and output fluxes
d_riv = 2 # riverine 
d_ht = 8 # hydrothermal vent
d_sr = 15 # subduction reflux 

# initialize values for Li sink into clays
d_sc0 <- 15
f_sc0 <- 2.9


# initialize values for Li isotopic balance in seawater
delSW0 <- 31

# initialize values for Li reservoir in seawater
M0 <- 3.4e6

# Li sink into clay differential for calculating isotope ratio of Li sink into clays (value used in Excel model)
sc_diff <- 16



#########      Define functions       #########

# Seawater Li Reservoir
fM <- function(dt.=dt, f_riv, f_ht, f_sr, f_sc) # t is timestep, M is seawater Li reservoir, and delSW is seawater Li isotopic balance
{
  dt*(f_riv + f_ht + f_sr - round(f_sc, 1))
}

# Li Isotopic Balance in Seawater
fdelSW <- function(dt.=dt, M, delSW, f_riv, d_riv, f_ht, d_ht, f_sr, d_sr, f_sc, d_sc) 
{
  dt/M*(f_riv*(d_riv - delSW) + f_ht*(d_ht - delSW) + f_sr*(d_sr - delSW) - f_sc*(d_sc - delSW))
}

# Li sink into clays flux; clay formation during seafloor alteration and sediment diagenesis
f_fsc <- function(f_sc00=f_sc0, M00=M0, dt.=dt, input_Li, output_Li){
  f_sc00*((input_Li - output_Li)*dt + M00)/M00
}

# isotope ratio of output (Li sink into clays) flux
f_dsc <- function(delSW, sc_diff){
  delSW - sc_diff
}


# initialize data frame to hold data
df <- list()

for (k in 1:length(friv)){
  
  df[[k]] <- data.frame(t=t)
  
  # add to data frame
  df[[k]]$f_riv <- friv[k]
  df[[k]]$f_ht <- f_ht
  df[[k]]$f_sr <- f_sr
  
  # add to data frame
  df[[k]]$d_riv <- d_riv
  df[[k]]$d_ht <- d_ht
  df[[k]]$d_sr <- d_sr
  
  df[[k]]$d_sc <- NA
  df[[k]]$f_sc <- NA
  df[[k]]$d_sc[1] <- d_sc0
  df[[k]]$f_sc[1] <- f_sc0
  
  df[[k]]$delSW <- NA
  df[[k]]$delSW[1] <- delSW0
  
  df[[k]]$M <- NA
  df[[k]]$M[1] <- M0
  
  # initialize values for cumulative input and output fluxes
  df[[k]]$input_Li <- NA
  df[[k]]$output_Li <- NA
  df[[k]]$input_Li[1] <- friv[k] + f_ht + f_sr
  df[[k]]$output_Li[1] <- f_sc0
  
}



###########    4th order Runge-Kutta Solver ###########


for (j in 1:length(friv)){
  
  for (i in 1:N_step)
  {
    # update time dependent variables
    #t_i <- df[[j]]$t[i]
    delSW_i <- df[[j]]$delSW[i]
    f_riv_i <- df[[j]]$f_riv[i]
    f_ht_i <- df[[j]]$f_ht[i]
    f_sr_i <- df[[j]]$f_sr[i]
    d_riv_i <- df[[j]]$d_riv[i]
    d_ht_i <- df[[j]]$d_ht[i]
    d_sr_i <- df[[j]]$d_sr[i]
    d_sc_i <- df[[j]]$d_sc[i]
    f_sc_i <- df[[j]]$f_sc[i]
    M_i <- df[[j]]$M[i]
    
    # calculate slopes
    m1 <- fM(dt, f_riv_i, f_ht_i, f_sr_i, f_sc_i)
    m2 <- fM(dt, f_riv_i, f_ht_i, f_sr_i, f_sc_i)
    m3 <- fM(dt, f_riv_i, f_ht_i, f_sr_i, f_sc_i)
    m4 <- fM(dt, f_riv_i, f_ht_i, f_sr_i, f_sc_i)
    
    k1 <- fdelSW(dt, M_i, delSW_i, f_riv_i, d_riv_i, f_ht_i, d_ht_i, f_sr_i, d_sr_i, f_sc_i, d_sc_i)
    k2 <- fdelSW(dt, M_i+m1/2, delSW_i+k1/2, f_riv_i, d_riv_i, f_ht_i, d_ht_i, f_sr_i, d_sr_i, f_sc_i, d_sc_i)
    k3 <- fdelSW(dt, M_i+m2/2, delSW_i+k2/2, f_riv_i, d_riv_i, f_ht_i, d_ht_i, f_sr_i, d_sr_i, f_sc_i, d_sc_i)
    k4 <- fdelSW(dt, M_i+m3, delSW_i+k3, f_riv_i, d_riv_i, f_ht_i, d_ht_i, f_sr_i, d_sr_i, f_sc_i, d_sc_i)
    
    # update dependent values
    M_new <- M_i + 1/6*(m1 + 2*m2 + 2*m3 + m4)
    delSW_new <- delSW_i+1/6*(k1 + 2*k2 + 2*k3 + k4)
    dsc_new <- f_dsc(delSW_new, sc_diff)
    fsc_new <- f_fsc(input_Li=df[[j]]$input_Li[i], output_Li=df[[j]]$output_Li[i])
    
    # update data frame
    df[[j]]$delSW[i+1] <- delSW_new
    df[[j]]$d_sc[i+1] <- dsc_new
    df[[j]]$f_sc[i+1] <- fsc_new
    df[[j]]$M[i+1] <- M_new
    df[[j]]$input_Li[i+1] <- df[[j]]$input_Li[i] + f_riv_i + f_ht_i + f_sr_i
    df[[j]]$output_Li[i+1] <- df[[j]]$output_Li[i] + fsc_new
    
  }
}


##########  plot numerical results ############

# remove scientific notation
options(scipen=999)

# plot delSW
plot(x=df[[1]]$t, y=df[[1]]$delSW, type="l", lwd=2, col="blue", main="d7Li (Seawater)", xlab="time", ylab="d7Li", ylim=c(10, 35))
lines(x=df[[2]]$t, y=df[[2]]$delSW, type="l", lwd=2, col="red")
lines(x=df[[3]]$t, y=df[[3]]$delSW, type="l", lwd=2, col="green")
#lines(x=df[[4]]$t, y=df[[4]]$delSW, type="l", lwd=2, col="orange")

# plot flux sink
plot(x=df[[1]]$t, y=df[[1]]$f_sc, type="l", lwd=2, col="blue", main="Flux Sink", xlab="time", ylab="flux", ylim=c(0, 25))
lines(x=df[[2]]$t, y=df[[2]]$f_sc, type="l", lwd=2, col="red")
lines(x=df[[3]]$t, y=df[[3]]$f_sc, type="l", lwd=2, col="green")
#lines(x=df[[4]]$t, y=df[[4]]$f_sc, type="l", lwd=2, col="orange")

