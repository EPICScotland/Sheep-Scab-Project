##############################################################################
#   Program for solving a differential equation based compartmental model for
#   sheep-scab.
##############################################################################

library('deSolve')
#library('rootSolve')
#library('bvpSolve')


##### Parameters

#\frac{dS}{dt}=-\frac{\beta}{N-R} S ((1-\phi) I_{sc}+\phi I_c) - (1-S_p) \pi(t) S +\delta P  $$
  
  #\frac{dI_{sc}}{dt}=\frac{\beta}{N-R} S ((1-\phi) I_sc+\phi I_c)-(1-\psi) I_{sc}+\rho I_c -S_e \pi(t) I_{sc}$$
  
  #\frac{dI_c}{dt}=(1-\psi) I_{sc}-\rho I_c-\mu I_c$$
  
  # \frac{dR}{dt}= \mu I_c$$

  #$$ \frac{dP}{dt}=  (1-S_p) \pi(t) S +S_e \pi(t) I_{sc} - \delta P $$

scab<- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    
    #Farm 1
    
    dS1 <-  - (beta / (N1-R1)) * S1 * ((1-gamma1)*((1 - phi) * I_sc1 + phi * I_c1) + gamma1*((1-phi)*I_sc2 +phi*I_c2)) - (1 - S_p) * T1 * S1 + delta * P1 # susceptibles
    dI_sc1 <- (beta / (N1-R1)) * S1 * ((1-gamma1)*((1 - phi) * I_sc1 + phi * I_c1) + gamma1*((1-phi)*I_sc2 +phi*I_c2))-(1 - psi) * I_sc1 + rho * I_c1 - S_e * T1 * I_sc1 #sub-clinicals
    dI_c1 <- (1 - psi) * I_sc1 - rho * I_c1 - mu * I_c1 #clinicals
    dP1 <-  (1 - S_p) * T1 * S1 + S_e * T1 * I_sc1 - delta * P1 #protecteds
    
    #Farm 2
    
    dS2 <-  - (beta / (N2-R2)) * S2 * ((1-gamma2)*((1 - phi) * I_sc2 + phi * I_c2) + gamma2*((1-phi)*I_sc1 +phi*I_c1)) - (1 - S_p) * T2 * S2 + delta * P2 # susceptibles
    dI_sc2 <- (beta / (N2-R2)) * S2 * ((1-gamma2)*((1 - phi) * I_sc2 + phi * I_c2) + gamma2*((1-phi)*I_sc1 +phi*I_c1)) - (1 - psi) * I_sc2 + rho * I_c2 - S_e * T2 * I_sc2 #sub-clinicals
    dI_c2 <- (1 - psi) * I_sc2 - rho * I_c2 - mu * I_c2 #clinicals
    dP2 <-  (1 - S_p) * T2 * S2 + S_e * T2 * I_sc2 - delta * P2 #protecteds
    
    #Combined removeds
    
    dR1 <- mu * I_c1  # removeds
    dR2 <- mu * I_c2
    
    res <- c(dS1, dI_sc1, dI_c1, dP1, dS2, dI_sc2, dI_c2, dP2, dR1, dR2)
    list(res)
  })
}
## The parameters
N01 <- 500
N02 <-500
days <- 365
beta0 <- N01*0.001
mu0 <- 0.0512933/days #following Jamie's suggestion for re-scaling mortality rate
rho0 <- 1/1000
parms <- c(N1 = N01, N2 = N02, beta = beta0, phi = 0.16875,  S_p = 0.965, S_e = 0.982, T1 = 0, T2 =0, delta = 0.1, psi = 1/40, 
           rho = rho0, mu = mu0, gamma1=0.0, gamma2 = 0.0)
## vector of timesteps
times <- seq(0, 365, by = 1) #You were outputing data 365 times a day!

## Start values for steady state
xstart <- c(S1 = N01, I_sc1 = 0, I_c1 = 0, P1 = 0, S2 = N02-1, I_sc2 = 1, I_c2 = 0, P2 = 0, R1 = 0, R2 = 0)
## Solve model

ptm<-proc.time()

out <- ode(y = xstart, times = times,
           func = scab, parms = parms) #solve ode system

duration<-proc.time()-ptm

## Default plot method

plot(out[,"time"],out[,"S1"],xlab="Time",ylab="S1,I_sc1,I_c1",pch=".", col="green",ylim=c(0,N01))
par(new=T)
plot(out[,"time"],out[,"I_sc1"],xlab="",ylab="",pch=".",yaxt="n", col="orange",ylim=c(0,N01))
par(new=T)
plot(out[,"time"],out[,"I_c1"],xlab="",ylab="",pch=".",yaxt="n", col="red",ylim=c(0,N01))
par(new=T)
plot(out[,"time"],out[,"P1"],xlab="",ylab="",pch=".",yaxt="n", col="black",ylim=c(0,N01))

plot(out[,"time"],out[,"S2"],xlab="Time",ylab="S2,I_sc2,I_c2",pch=".",col="green",ylim=c(0,N02))
par(new=T)
plot(out[,"time"],out[,"I_sc2"],xlab="",ylab="",pch=".",yaxt="n", col="orange",ylim=c(0,N02))
par(new=T)
plot(out[,"time"],out[,"I_c2"],xlab="",ylab="",pch=".",yaxt="n", col="red",ylim=c(0,N02))
par(new=T)
plot(out[,"time"],out[,"P2"],xlab="",ylab="",pch=".",yaxt="n", col="black",ylim=c(0,N02))
#legend("bottomright")

plot(out[,"time"],out[,"S2"],xlab="Time",ylab="S2,I_sc2,I_c2",type="l",col="green",ylim=c(0,N02))
lines(out[,"time"],out[,"I_sc2"],xlab="",ylab="",pch=".",yaxt="n", col="orange")
points(out[,"time"],out[,"I_c2"],xlab="",ylab="",yaxt="n", col="red")
lines(out[,"time"],out[,"P2"],xlab="",ylab="",pch=".",yaxt="n", col="black")

#plot(out)

