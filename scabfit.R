#######################################
#
# Program to estimate parameters of sheep-scab for the UK
# The idead is to try out the inverse modelling package FME following
# Soetart's paper on inverse modelling.
# In addition will try Grenfelll's underreporting estimation method in conjunction with this.
#
#######################################
library('deSolve')
library('FME')
library('xts')

# Prep the data
rawdat <-read.csv("scabData.csv",header=TRUE)
dat <- data.frame(cbind(rawdat[,8],rawdat[,7],rawdat[,6],rawdat[,5],rawdat[,4],rawdat[,3]))
ts <- stack(dat)
C <- ts[,1]
time <-c(1:72)
C <- cbind(time,C)
cases <- C[,2]
ts.plot(cases) # plot time seris

# Note following Finkenstadt and Grenfell 2000 I=theta C (C cases)
# Assume only clinical cases reported the Ic = c * C and Isc = C- c * C

N0    <- 100
days <- 365
years <- 1
Timehorizon  <- days*years
# delta <- 0.01 
beta0 <- N0*0.0083 # where 0.0083? Re-scaled 
phi0  <- 0.16875
gamma <- 0.0 # 0.001
Sp   <- 0.965
Se <- 0.982
mu0 <- 0.0512933/days #following Jamie's suggestion for re-scaling mortality rate
rho0 <- 1-mu0
N <- N0
c <- 1 # proportion of observed cases
theta <- 2*pi/(Timehorizon*0.4) # period



scab <- function(t, y, parms) {
  with(as.list(c(parms, y)), {
    
    C <- Ic/c
      
    dS <- (- S * beta * cos(theta*Z)* (Isc + c * C * phi)/(N - R) ) 
    
    dIsc <-  (c * C  * rho - Isc * psi  + S * beta * cos(theta*Z) * (Isc + c * C * phi)/(N - R) )
    
    dIc <- (Isc * psi - c * C * mu - c * C * rho ) 
    
    dR <-  c * C * mu 
    
    dC <- (Isc * psi - c * C * mu - c * C * rho )/c
    
    dZ <- 1
    
    
    res <- c(dS, dIsc, dIc, dR, dC, dZ)
    list(res)
  })
}

# solve odesystem

parms <- c(N = N0, beta = beta0, phi = 0.16875/days, psi = 1/40, rho = rho0, mu = mu0, c=c,theta=theta)
## vector of timesteps
times <- seq(0,Timehorizon, by = 0.1) # 0,5 step gives good results

xinit <- c(S = N0-1, Isc = 1, Ic = 0, R = 0, C= ts[1,1], Z=0)
out <- ode(y = xinit, times = times,
            func = scab, parms = parms,outnames=NULL)

# residual calculations
# parms <- c(N = N0, beta = beta0, phi = 0.16875/days, psi = 1/40, rho = rho0, mu = mu0, c=c,theta=theta)

scabcost <- function (P) {
  parms["N"] <-  P[1]
  parms["beta"] <- P[2]
  parms["phi"] <-  P[3]
  parms["psi"] <- P[4]
  parms["rho"] <- P[5]
  parms["mu"] <-  P[6]
  parms["c"] <- P[7]
  parms["theta"] <- P[8]
  out <-  ode(y = xinit, times = times,
              func = scab, parms = parms,outnames=NULL,maxsteps=8000)
  #cost <- modCost(scab,ts[,1],x="time",error="sd")
  return(modCost(out,obs=C,error="sd"))
}

# Sensitivity/Identifiability analysis

#Coll <- collin(sF <- sensFun(func = scabcost, parms = parms, varscale = 1))
#plot(Coll, log = "y")
#abline(h = 20, col = "red")
#collin(sF,parset=2:3)

# I need to work on identification of the model parameters a buit more.
# Fit the model

MCMC <- modMCMC(f=scabcost,p=parms,lower=c(0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0))

#MCMC <- modMCMC(f = scabcost, p = Fit$par, niter = 5000,
 #                wvar0 = 0.1, updatecov = 50)

Fit <- modFit(f = scabcost, p = parms, 
              lower=c(0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0),method="Nelder-Mead") #This worked

# Produces a fit but not great the data is very noisy

summary(Fit)
Status API Training Shop Blog About Pricing
© 2016 GitHub, Inc. Terms Privacy Security Contact Help
