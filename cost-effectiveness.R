#Cost effectiveness analysis of sheep scab using markov chains
#Have added DALY's as the outcome variable also discounting of costs based on each iteration representing a day
# Timescale is calibrated against incubation period of scab from subclinical until clinical signs appear
#this is based on experimental studies. I have calculated DALY's cumulatively

rm(list=ls())
library('markovchain')
library('igraph')
library('BCEA')
library('MASS')

#some notes on mite lifecycle: females 1:3 eggs/day adults live 50 days and 
#lay 50-100 eggs. Egg to egg lifecycle 10-14 days minimum. Can survive 2-3weeks off hostup to 12 weeks in winter.
# subclinical period lasts about 10-35 days (experimental) Need this to calibrate discount rate.

#Parameters
replications<-1 #Number of replications of monte-carlo and CEA to do.
Years<-1.1
iterations<-365*Years
numtreatments<-1
testcost<-12 #8 from AHVLA
scrapingcost<-15 #Easter Ross Vets but confirmed by SRUC consulting £15 per scraping 4-8 scrapings needed
r<-0.03
delta<-1/(1+r/365) #discount factor daily
k<-25000 #wtp
weightgain<-0.1317 #kg per day
birthweight<-4.4
treatmentdate<-(iterations-70)/numtreatments
flocksize<-4

#The mean birth-weight for lambs in dataset 1 was 4.4±1.0 kg. The mean birth-weight for each farm 
#ranged from 3.8±0.8 to 5.1±1.0 kg. The mean birth-weight for lambs in dataset 2 was 4.6±1.1 kg. 
#The mean birth-weight for each farm ranged from 4.0±0.9 to 5.1±1.0 kg.


#Sensitivity and Specificity

se<-0.93
sp<-0.9

#Cost parameters

treatmentcost<-1 #from NAAC tables
weightloss<-0.4*0.53+0.49*0.1 #0.4 proportion of ewes 2012 0.49 proportion of lams 2012 Weightloss for rams not available.


#the following just generates and runs a Markov chain of sheep-scab ignoring treatment

#Set-up Markov chain

tune<-0.95
mort<-0.05

P<-matrix(c(1-0.0083,0.0083,0.0,0.0,
            0.0,1-0.16875,0.16875,0.0,
            0.0,0.95,0.0,0.05, #choose probability of remaining clinically so that 60-90 days is the clinical ok but this indicates mortality rates is too high duration
            0.0,0.0,0.0,1.0),byrow=T,nrow=4)

diseaseStates<-c("S","I_sc","I_c","D")

mcDisease<-new("markovchain",states=diseaseStates,byrow=TRUE,transitionMatrix=P,name="scab")

Igraph<-as(mcDisease, "igraph")

#repeat loop
result=list()
for (i in 1:replications){
#simulates the chain
chain<-rmarkovchain(n = iterations, object = mcDisease, t0 = "S") 

#Calculate probability of of didease on day 4
firstPassagePdF <- firstPassage(object = mcDisease, state = "I_sc", n = iterations)
hit<-firstPassagePdF[4,1]

#print some stuff

print(P)
print(chain)

#creates a frequency table of the simulated chain transitions

freq<-createSequenceMatrix(chain)
print(freq)

#interactive plot of the state transitions trying to figure oput how to label the edges 
#with the transition probabilities

#mcIgraph<-as(mcDisease, "igraph")
plot(mcDisease)

#Incorporating a 50:50 chance of treatment in sub-clinically infected and susceptible states
#Assume treatment succesful in clinically infected (would probably cull here) treament would mean transition with prob 1 to
# S cull would mean transition with prob 1 to D

#transitions probs given treating

PT<-matrix(c(1,0,0,0,
             1,0,0,0,
             1,0,0,0,
             0,0,0,1),byrow=T,nrow=4)

#transition probs 50:50 chance of treatment 

s<-matrix(c(0,0,0,0, #first entry 0 not 0.5 because treatment can't impact a healthy flock
            0,0.5,0,0,
            0,0,1,0,
            0,0,0,1),byrow=T,nrow=4)

PTNotest<-s%*%PT+(diag(4)-s)%*%P

mcDiseaseNoTest<-new("markovchain",states=diseaseStates,byrow=TRUE,
                   transitionMatrix=PTNotest,name="scab")

chainNoTest<-rmarkovchain(n = iterations, object = mcDiseaseNoTest, t0 = "S")


#matrix contains treatment probabilities in different states given test has been adopted

#st<-matrix(c(1-sp,0,0,0,
 #             0,1-se,0,0,
  #             0,0,1,0,
   #            0,0,0,1),byrow=T,nrow=4)

st<-matrix(c(0,0,0,0,
             0,se,0,0,
             0,0,1,0,
             0,0,0,1),byrow=T,nrow=4)

PTtest<-st%*%PT+(diag(4)-st)%*%P

#transitions after treatment during rotection period

PP<-matrix(c(1,0,0,0,
             1,0,0,0,
             1,0,1,0,
             0,0,0,1),byrow=T,nrow=4)

mcDiseaseTest<-new("markovchain",states=diseaseStates,byrow=TRUE,transitionMatrix=PTtest,name="scab")

Igraph<-as(mcDiseaseTest, "igraph")

#simulates the chain
chainTest<-rmarkovchain(n = iterations, object = mcDiseaseTest, t0 = "S") 

print(chainTest)

#configuration

statesTest<-as.vector(chainTest)
statesNoTest<-as.vector(chainNoTest)
costT<-matrix(NA,nrow=iterations,ncol=2)
outcome<-matrix(NA,nrow=iterations,ncol=2)
weightT<-matrix(NA,nrow=iterations,ncol=1)
weightNT<-matrix(NA,nrow=iterations,ncol=1)
rownames(costT)<-c(1:iterations)
colnames(costT)<-c("T","NT")
rownames(outcome)<-c(1:iterations)
colnames(outcome)<-c("T","NT")
incidencedT<-0 #counter for disability events with test
incidencedNT<-0 #counter for disability events without test
incidencelT<-0 #counter for duration of death with test 
incidencelNT<-0 #counter for duration of death with no test 
weight<-birthweight

#populate test case

for (state in 1:iterations) {
  t<-state-1
  weight<-weight+weightgain
  
  if (statesTest[state]=="S") {
    weightT[state,1]<-weight  
  if (t==treatmentdate)  
  costT[state,1]<-delta^t*(testcost+(1-sp)*treatmentcost)
  else
  costT[state,1]<-0.0  
  #outcome[state,1]<-1
  #costT[state,2]<-(1-sp)*treatmentcost
 
  }
  else 
  {
  if (statesTest[state]=="I_sc") {
    costT[state,1]<-delta^t*(testcost+se*treatmentcost) 
    #costT[state,2]<-se*treatmentcost
    #outcome[state,1]<-weightloss
    incidencedT<-incidencedT+1
    weightT[state,1]<-weight-weight*0.53/incidencedT
  }
  else
  {
  if (statesTest[state]=="I_c") {
    costT[state,1]<-delta^t*treatmentcost
    #costT[state,2]<-treatmentcost
    #outcome[state,1]<-weightloss
    incidencedT<-incidencedT+1
    weightT[state,1]<-weight-weight*0.53/incidencedT
  }
  else 
    costT[state,1]<-0
    #costT[state,2]<-0
    #outcome[state,1]<-0
   incidencelT<-incidencelT+1
   weightT[state,1]<-0.0
  
  }
}
yldT<-incidencedT*weightT[state,1]
yllT<-incidencelT  
outcome[state,1]<-(yldT+yllT)/365 #(DALY version)
#need to apportion weight loss over duration of illness
}                     

weight<-birthweight
for (state in 1:iterations) {
  weight<-weight+weightgain
  if (statesNoTest[state]=="S") {
    #costT[state,1]<-testcost+(1-sp)*treatmentcost
    costT[state,2]<-delta^t*(scrapingcost+0.5*treatmentcost)
    #outcome[state,2]<-1
    weightNT[state,1]<-weight
  }
  else 
  {
    if (statesNoTest[state]=="I_sc") {
      #costT[state,1]<-testcost+se*treatmentcost 
      costT[state,2]<-delta^t*(scrapingcost+0.5*treatmentcost)
      #outcome[state,2]<-weightloss
      incidencedNT<-incidencedNT+1
      weightNT[state,1]<-weight-weight*0.53/incidencedNT
    }
    else
    {
      if (statesNoTest[state]=="I_c") {
        #costT[state,1]<-treatmentcost
        costT[state,2]<-delta^t*treatmentcost  
        #outcome[state,2]<-weightloss
        incidencedNT<-incidencedNT+1
        weightNT[state,1]<-weight-weight*0.53/incidencedNT
      }
      else 
      #costT[state,1]<-0
      costT[state,2]<-0
      #outcome[state,2]<-0
     incidencelNT<-incidencelNT+1
     weightNT[state,1]<-0.0
      
    }
  }
  yldNT<-incidencedNT*weightNT[state,1]
  yllNT<-incidencelNT  
  outcome[state,2]<-(yldNT+yllNT)/365
}          


#Outcomes


ints<-c("Adopt test","Don't adopt test")
result[[i]]<-bcea(outcome, costT, interventions = ints, Kmax = 25000, plot = TRUE)
#summary(result[[i]],wtp=k) #am using the approximate value of a lamb as the WTP measure current prices £167, ewes 56.79
# if we take a  proportional weight based on population proportions to approximate societal value 0.4*56.79+0.49*167=104 and round to 100
# wtp must be increments of 100 then this seems reasonable
summary(result[[i]])
#plot(result[[i]],wtp=k)
plot(result[[i]])
}
#contour2(result)
#ceplane.plot(result)
#ceplane.plot(result, # plots the Cost-Effectiveness plane
 #            comparison=1, # if more than 2 interventions, selects the
             # pairwise comparison
  #           wtp=25000, # selects the relevant willingness to pay
             # (default: 25,000)
   #          graph="base" # selects base graphics (default)
#)
