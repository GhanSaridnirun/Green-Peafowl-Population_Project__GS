library(nimble)

mySeed <- 0
set.seed(mySeed)


# Data Bundle

#Month order and season for data

# month <- factor(c("Nov19","Dec19",
#                   "Jan20","Feb20","Mar20","Apr20","May20","Jun20","Jul20","Aug20",
#                   "Sep20","Oct20","Nov20","Dec20","Jan21","Feb21","Mar21","Apr21",
#                   "May21","Jun21","Jul21","Aug21","Sep21","Oct21","Nov21","Dec21"
# ),ordered = TRUE)
# 
# season <- factor(c("Breeding","Breeding","Breeding","Breeding","Breeding",
#                    "Breeding","Non","Non","Non","Non","Non","Non","Breeding","Breeding",
#                    "Breeding","Breeding","Breeding","Breeding","Non","Non","Non","Non",
#                    "Non","Non","Breeding","Breeding"
# ),ordered = TRUE)


# Count in Non-Breeding Season

ChF_NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2)             #Season Label for chicks count
ChF_NB <- c(16,26,69,80,60,20,25,24,48,59,11,12)    #Female Chicks Count

AF_NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2)              #Season Label for breeder count
#Br_NB <- c(35,34,83,78,56,25,29,28,47,57,14,8)      #Breeder {Female with chicks} count in Non-Breeding 
AF_NB <- c(71,69,134,108,86,33,101,68,71,90,34,13)  #All Female 


# Count in Breeding Season

JuF_BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)          #Season Label for juvenile count
JuF_BN <- c(9,18,33,33,10,13,6,9,13,27,14,25,11,0)   #Female Juvenile count

AF_BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)           #Season Label for breeder count
#Br_BN <- c(9,23,32,25,12,10,14,9,11,20,13,27,0,0)    #Breeder {Female with juveniles} count in Breeding
AF_BN <- c(18,49,79,71,47,44,30,25,38,54,48,97,17,1) #All Female



# # Single Female count data
# SF_BN <- c(9,26,48,46,35,34,16,16,27,34,35,70,9,1)   #Single Female count in Breeding 
# SF_NB <- c(36,35,51,30,30,8,72,40,24,33,20,5)        #Single Female count in Non-Breeding         
# 
# 
# # Male count Data
# 
# ChM_NB <- c(20,39,78,72,52,21,25,26,33,54,13,4)      #Male Chicks Count
# JuM_BN <- c(13,17,28,29,8,8,20,9,10,36,11,19,11,0)   #Male juvenile count
# 
# M1y_BN <- c(4,16,11,19,13,18,8,2,9,14,16,28,1,0)     #Male 1 year count in Breeding   
# M1y_NB <- c(9,6,9,13,5,5,12,19,14,10,5,10)           #Male 1 year count in Non-Breeding
# 
# M2y_BN <- c(1,6,4,6,2,6,4,6,1,8,2,14,0,0)            #Male 2 years count in Breeding
# M2y_NB <- c(5,1,8,4,0,2,2,1,0,0,0,1)                 #Male 2 years count in Non-Breeding
# 
# M3y_BN <- c(27,118,164,119,49,57,54,60,76,91,111,83) #Male 3 years count in Breeding
# M3y_NB <- c(75,57,143,139,85,34,71,63,83,82,56,37)   #Male 3 years count in Non-Breeding

ny.data <- 3 # Number of years for which the data collected
ny.sim <- 0 # Number of years to simulate after the data collection

## Arrange constants

GP.IPMconstants <- list(Tmax = ny.data + ny.sim, 
                        Amax = 4,
                        ny.sim = ny.sim,
                        ny.data = ny.data,
                        ChF_NB_yr = ChF_NB_yr,
                        AF_NB_yr = AF_NB_yr,
                        JuF_BN_yr = JuF_BN_yr,
                        AF_BN_yr = AF_BN_yr
                        )

## Arrange data
GP.IPMdata <- list(ChF_NB = ChF_NB,  
                   AF_NB = AF_NB,  
                   JuF_BN = JuF_BN,  
                   AF_BN = AF_BN 
                   ) 


GPIPM <- list(GP.IPMconstants, GP.IPMdata)  
str(GPIPM)  
  
  
## Nimble code for run the whole model

GP.IPMcode <- nimbleCode({
  
  
  # Priors and linear models
  
  # Survival 
  
  s_NB[1] ~ dunif(0.30, 0.40)  # chick (half year)
  s_BN[1] ~ dunif(0.60, 0.80)  # juvenile (half year)
  
  s_yr_sa ~ dunif(0.50, 0.70) # sub adult (whole year)
  s_NB[2:3] <- sqrt(s_yr_sa) # sub adult (half year, N -> B)
  s_BN[2:3] <- sqrt(s_yr_sa) # sub adult (half year, B -> N)
  
  s_yr_ad ~ dunif(0.60, 0.80)  # breeder (whole year)
  s_NB[4] <- sqrt(s_yr_ad) # breeder (half year, N -> B)
  s_BN[4] <- sqrt(s_yr_ad) # breeder (half year, B -> N)
  
  
  
  # Productivity
  
  for (t in 1:ny.data){
    log.rho[t] ~ dnorm(log(mean.rho), sd = sigma.rho)
    rho[t] <- exp(log.rho[t])
  }
  
  mean.rho ~ dunif(3, 5)
  sigma.rho ~ dunif(0, 1) 
  
  # Sex ratio of the chicks
  
  gamma ~ dbeta(1, 1)  
  
  
  #Detection Probability
  
  # for (t in 1:ny.data){
  #   logit.p[t] ~ dnorm(lmean.p, sd = sigma.p)
  #   p[t] <- ilogit(logit.p[t])
  # }
  # 
  # mean.p ~ dbeta(1, 0.9)
  # lmean.p <- logit(mean.p)
  # sigma.p ~ dunif(0, 5)
  
  for (t in 1:ny.data){
    p[t] <- 0.9 # TODO update initial vales for p
  }
  
  # Population count data (state-space model)
  # Model for the initial population size in Breeding season
  
  for(a in 1:4){
    pinit[a] ~ dunif(0, 200)
    NBreed[a,1] ~ dpois(pinit[a]) 
  }
  
  # Book-keeping
  
  NNon[1:Amax,1] <- 0   
  
  
  
  # Process model: Breeding -> Non-Breeding season transition
  
  for (t in 1:ny.data){
    
    # Total number of chicks
    
    Fec[t] ~ dpois(sum(NBreed[3:4,t]) * rho[t])
    
    
    # Allocate chicks to a sex
    
    NNon[1,t+1] ~ dbin(gamma, Fec[t]) # Female chicks 
    
    
    # Survival
    
    for(a in 2:3){
      NNon[a,t+1] ~ dbin(s_BN[a-1], NBreed[a-1,t])
    }
    
    NNon[4,t+1] <- surv_NBreed3[t+1] + surv_NBreed4[t+1]
    surv_NBreed3[t+1] ~ dbin(s_BN[3], NBreed[3,t])
    surv_NBreed4[t+1] ~ dbin(s_BN[4], NBreed[4,t])
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:4){
      NBreed[a,t+1] ~ dbin(s_NB[a], NNon[a,t+1])
    }
      
    
    # Observation Model in Non-Breeding
    
      for(j in 12) {
      
      ChF_NB[j] ~ dpois(p[ChF_NB_yr[j]] * NNon[1,ChF_NB_yr[j]])
      AF_NB[j] ~ dpois(p[AF_NB_yr[j]] * sum(NNon[2:4,AF_NB_yr[j]]))
      
    }
    
    # Observation Model in Breeding
    
      for(h in 14) {
      
      JuF_BN[h] ~ dpois(p[JuF_BN_yr[h]] * NBreed[1,JuF_BN_yr[h]])
      AF_BN[h] ~ dpois(p[AF_BN_yr[h]] * sum(NBreed[2:4,AF_BN_yr[h]]))
      
      }
    } 
  }
)

# Initial values

source("GPeafowlIPM_InitialSim.R")

  Inits <- GP_IPM_Init(Tmax = ny.data + ny.sim)
Inits


# Parameters monitored
parameters <- c("s_NB", "s_BN", "mean.rho", "gamma", "sigma.rho", "rho", "p",
                "NBreed", "NNon", "Fec")


# MCMC settings

ni <- 10
nb <- 0
nt <- 1
nc <- 3

# ni <- 200000     # Run Time around 15 minutes
# nb <- 50000      # Using initial values
# nt <- 30
# nc <- 4

# ni <- 500000     # Run Time around 30 minutes
# nb <- 1000      # Using initial values
# nt <- 30
# nc <- 4

#ni <- 1000
#nb <- 0
#nt <- 30
#nc <- 4



# Run
out <- nimbleMCMC(code = GP.IPMcode,
                  data = GP.IPMdata,
                  constants = GP.IPMconstants,
                  inits = Inits,
                  monitors = parameters,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = nb,
                  niter = ni,
                  thin = nt,
                  nchains = nc)
  
print(out,2)

library(MCMCvis)

  MCMCsummary(out, params = 'NBreed', round = 2)

  MCMCsummary(out, params = 'NNon', round = 2)

  MCMCsummary(out, params = 'Fec', round = 2)

  MCMCsummary(out, params = 'rho', round = 2)

  MCMCsummary(out, params = 'p', round = 2)

  MCMCsummary(out, 
            params = c('s_NB','s_BN','gamma','mean.rho',
                       'sigma.rho','mean.p','sigma.p'),
            round = 2)

# Trace Plot
  
  MCMCtrace(out, params = 'NBreed',
          pdf = TRUE,
          Rhat = TRUE,
          n.eff = TRUE)


  MCMCtrace(out, params = 'NNon',
            pdf = TRUE,
            Rhat = TRUE,
            n.eff = TRUE)
  
  




