#Gillespie algorithm implementation of sheep scab
# dS/dt=-beta*S((1-phi)*I_sc+phi*I_c)
# dI_sc/dt = beta*(1-phi)*S*I_sc+beta*phi*I_c*S-phi*I_sc+rho*I_c
#dI_c/dt=phi*I_sc-rho*I_c-mu*I_c have introduced something here to remain in clinical class will see if it works
#dM/dt=mu*I_c
#protection equation
#dP/dt=pi(t)*S+pi*I_sc
#dw/dt=w-l*I_sc
#pi(t)
#phi is the proportion of clinicals
#beta0 background infection from AHVLA
#Trying tmake beta time dependent with seasonal forcing may need to add demography to get this to work.
#recovery rate rho=1-mu, mu is mortality rate
library("GillespieSSA")

Years<-1
days<-365*Years
N<-1000 #no. of sheep
beta0<-N*0.0083 #re-scale infection
mu0<-0.0512933/days #following Jamie's suggestion for re-scaling mortality rate
rho0<-1-mu0
theta<-2*pi*(2/days) #seasonal forcing parameter
parms<-c(beta=beta0,phi=0.16875,mu=mu0,rho=rho0,delta=0.5,sp=0.965,se=0.982,pi=1) #phi is clinical proortion delt is prtective decay parameter epsilon is immune proportion not being used. rho recovery parameter.
#x0<-c(S=95,I_sc=5,I_c=0,M=0,t=0)
x0<-c(S=999,I_sc=1,I_c=0,M=0,P=0)
tune<-(-1.1)
shift<-(-1.5)

#state update matrix
#nu<-matrix(c(-1,+1,-1,0,0,0,
 #            +1,-1,+1,-1,+1,0,
  #           0,0,0,+1,-1,-1,
   #          0,0,0,0,0,+1
    #          ),nrow=4,byrow=T)

#state update matrix with treatment

nu<-matrix(c(-1,+1,-1,-1,+1,+1,0,0,0,0,0, #S done
            +1,-1,+1,0,0,0,-1,+1,+1,-1,0, #I_sc check
           0,0,0,0,0,0,+1,-1,-1,0,-1,#I_c ok
           0,0,0,0,0,0,0,0,0,0,+1,#M ok
           0,0,0,+1,+1,-1,0,0,0,0,0 #P ok
          ),nrow=5,byrow=T)

a<-c("(beta/(N+1-M))*S*I_sc",
     "(beta/(N+1-M))*phi*S*I_sc",
     "(beta/(N+1-M))*phi*S*I_c",
     "pi*S",
     "pi*sp*S",
     "delta*P",
     "I_sc",
     "phi*I_sc",
     "rho*I_c",
     "se*pi*I_sc",
     "mu*I_c")

#state update matrix with seasonal forcing and demography

#nu<-matrix(c(+1,-1,+1,-1,0,0,0,0,
 #            0,+1,-1,+1,-1,+1,0,0,
  #           0,0,0,0,+1,-1,-1,0,
   #          0,0,0,0,0,0,+1,0,
    #         0,0,0,0,0,0,0,+1
#),nrow=5,byrow=T)


#a<-c("(beta/(N+1-M))*S*I_sc","(beta/(N+1-M))*phi*S*I_sc","(beta/(N+1-M))*phi*S*I_c","phi*I_sc",
#     "rho*I_c","mu*I_c") #Forces of infection etc.

#with treatment
#a<-c("(beta/(N+1-M))*S*I_sc","(beta/(N+1-M))*phi*S*I_sc","(beta/(N+1-M))*phi*S*I_c","epsilon*S","pi*S","pi*sp*S","delta*P",
 #    "I_sc","phi*I_sc","rho*I_c","se*pi*I_sc","mu*I_c","A_s","A_sc")



#rates with seasonal forcing
#a<-c("mu*(1+tune*sin(theta*t+shift))*N","((beta*(1+tune*sin(theta*t)))/(N+1-M))*S*I_sc",
 #    "((beta*(1+tune*sin(theta*t)))/(N+1-M))*phi*S*I_sc",
  #   "((beta*(1+tune*sin(theta*t)))/(N+1-M))*phi*S*I_c","phi*I_sc",
   #  "rho*I_c",
    # "mu*I_c",
     #"1")

out <- ssa(x0, a, nu, parms, tf = days, method="D",
            simName = "sheep-scab",censusInterval=1)

#out2<-ssa.run(x0, a, nu, parms,tau=1,f=1.5,epsilon=0.01,nc=1,tf = 365, method="D")

ssa.plot(out)

testdate<-(days/4)*3/2 #treament set to middle of third quarter of year day 136
#note 70 day meat withdrawal means treating no later than day 299 robably earlier depending on animal weight
dayafter<-testdate+1
#pre-test calculations
P<-matrix(NA,nrow=days+1,ncol=1)
S<-matrix(NA,nrow=days+1,ncol=1)
I_sc<-matrix(NA,nrow=days+1,ncol=1)
P[1:testdate-1,1]<-out$data[1:testdate-1,6]-(1-sp)*pi*out$data[1:testdate-1,2]-se*out$data[1:testdate-1,2]
S[1:testdate-1,1]<-out$data[1:testdate-1,2]+(1-sp)*pi*out$data[1:testdate-1,2]
I_sc[1:testdate-1,1]<-out$data[1:testdate-1,3]+se*out$data[1:testdate-1,2]
#test date calculations
P[testdate,1]<-out$data[testdate,6]
S[testdate,1]<-out$data[testdate,2]
I_sc[testdate,1]<- out$data[testdate,3]
#post-test calculations
P[dayafter:days+1,1]<-out$data[dayafter:days+1,6]-(1-sp)*pi*out$data[dayafter:days+1,2]-se*out$data[dayafter:days+1,2]
S[dayafter:days+1,1]<-out$data[dayafter:days+1,2]+(1-sp)*pi*out$data[dayafter:days+1,2]
I_sc[dayafter:days+1,1]<- out$data[dayafter:days+1,3]+se*out$data[dayafter:days+1,2]




#Treament at middle of thid quarterso (days/4)*3/2

#plot(out$data[,6],out$data[,3],xlab="time",ylab="No.",cex=0.5,col="yellow")
#par(new=T)
#plot(out$data[,6],out$data[,4],xlab="",ylab="",col="red",cex=0.5,yaxt="n")
#par(new=T)
#plot(out$data[,6],out$data[,5],xlab="",ylab="",col="black",cex=0.5,yaxt="n")
#par(new=T)
#plot(out$data[,6],out$data[,2],xlab="",ylab="",col="green",cex=0.5,yaxt="n")

#steps<-dim(out$data)

#numsteps<-steps[1]

#timesteps<-matrix(NA,nrow=numsteps,ncol=1)
#for (i in 1:numsteps){
#timesteps[i]<-(out$data[i,1]-out$data[i-1,1])/numsteps
#}
Status API Training Shop Blog About Pricing
© 2016 GitHub, Inc. Terms Privacy Security Contact Help
