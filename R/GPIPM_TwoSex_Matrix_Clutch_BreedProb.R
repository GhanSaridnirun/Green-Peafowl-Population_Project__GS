
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

# Year Label
NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2) + 1 

#ChF_NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2) + 1        #Season Label for chicks count
ChF_NB <- c(16,26,69,80,60,20,25,24,48,59,11,12)    #Female Chicks Count

#AF_NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2) + 1              #Season Label for breeder count
#Br_NB <- c(35,34,83,78,56,25,29,28,47,57,14,8)      #Breeder {Female with chicks} count in Non-Breeding 
AF_NB <- c(71,69,134,108,86,33,101,68,71,90,34,13)  #All Female 

# Count in Breeding Season

# Year Label

BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)

#JuF_BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)          #Season Label for juvenile count
JuF_BN <- c(9,18,33,33,10,13,6,9,13,27,14,25,11,0)   #Female Juvenile count

#AF_BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)           #Season Label for breeder count
#Br_BN <- c(9,23,32,25,12,10,14,9,11,20,13,27,0,0)    #Breeder {Female with juveniles} count in Breeding
AF_BN <- c(18,49,79,71,47,44,30,25,38,54,48,97,17,1) #All Female

# # Single Female count data
# SF_BN <- c(9,26,48,46,35,34,16,16,27,34,35,70,9,1)   #Single Female count in Breeding 
# SF_NB <- c(36,35,51,30,30,8,72,40,24,33,20,5)        #Single Female count in Non-Breeding         

# Brood Data

# 1ch <- c(3,15,19,4,8,7,19,12,47,39,34,18,8,6,6,8,8,17,11,15,27,25,11,3,1,0)
# 2ch <- c(2,9,13,13,3,3,15,22,33,26,15,7,4,3,5,7,4,8,15,10,14,20,2,2,3,0)
# 3ch <- c(2,2,2,7,0,0,8,5,13,16,7,2,2,1,3,11,1,3,4,3,5,13,2,3,2,0)
# 4ch <- c(2,2,5,3,1,1,2,2,5,5,7,1,0,1,1,4,1,2,5,1,2,3,0,0,2,0)
# 5ch <- c(0,0,0,1,0,1,1,0,0,1,2,1,1,0,0,2,1,0,0,1,2,1,1,0,0,0)
# 6ch <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0)
# 7ch <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)
# brood <- c(9,28,39,28,12,12,45,41,98,87,65,29,15,11,15,33,15,30,36,31,50,62,16,8,8,0)

# Clutch Sizes

# ClutchSize <- c(3,7,4,4,4,11,3,3,5,5,10)

# # Male count Data
# 
ChM_NB <- c(20,39,78,72,52,21,25,26,33,54,13,4)      #Male Chicks Count
JuM_BN <- c(13,17,28,29,8,8,20,9,10,36,11,19,11,0)   #Male juvenile count

M1y_BN <- c(4,16,11,19,13,18,8,2,9,14,16,28,1,0)     #Male 1 year count in Breeding
M1y_NB <- c(9,6,9,13,5,5,12,19,14,10,5,10)           #Male 1 year count in Non-Breeding

M2y_BN <- c(1,6,4,6,2,6,4,6,1,8,2,14,0,0)            #Male 2 years count in Breeding
M2y_NB <- c(5,1,8,4,0,2,2,1,0,0,0,1)                 #Male 2 years count in Non-Breeding

M3y_BN <- c(27,118,164,119,49,57,54,60,76,91,111,83,26,1) #Male 3 years count in Breeding
M3y_NB <- c(75,57,143,139,85,34,71,63,83,82,56,37)   #Male 3 years count in Non-Breeding

M_BN <- rbind(JuM_BN, M1y_BN, M2y_BN, M3y_BN)
M_NB <- rbind(ChM_NB, M1y_NB, M2y_NB, M3y_NB)

ny.data <- 3 # Number of years for which the data collected
ny.sim <- 5 # Number of years to simulate after the data collection

# Reproduction data

source("R/ReproductionDataPrep.R")

## Arrange constants

GP.IPMconstants <- list(Tmax = ny.data + ny.sim, 
                        Amax = 4,
                        ny.sim = ny.sim,
                        ny.data = ny.data,
                        NB_yr = NB_yr,
                        BN_yr = BN_yr,
                        xmax = xmax,
                        ymax = ymax,
                        Year_BS = Year_BS 
                        
)

## Arrange data
GP.IPMdata <- list(ChF_NB = ChF_NB,
                   JuF_BN = JuF_BN,
                   AF_NB = AF_NB,
                   AF_BN = AF_BN,
                   M_BN = M_BN,
                   M_NB = M_NB,
                   Rep = Rep,
                   BroodSize = BroodSize
) 


GPIPM <- list(GP.IPMconstants, GP.IPMdata)  
str(GPIPM)  





## Nimble code for run the whole model

GP.IPMcode <- nimbleCode({
  
  
  # Priors and linear models
  
  # Survival 
  
  ## Chicks and juveniles
  
  sF_NB[1] ~ dunif(0.30, 0.40) 
  sF_BN[1] ~ dunif(0.60, 0.80) 
  
  sM_NB[1] <- sF_NB[1] 
  sM_BN[1] <- sF_BN[1]  
  
  ## Adults (no sex difference)
  
  s_yr_sa ~ dunif(0.50, 0.70) 
  s_yr_ad ~ dunif(0.60, 0.80)  
  
  s_yr_saF <- s_yr_sa
  s_yr_saM <- s_yr_sa 
  
  s_yr_adF <- s_yr_ad
  s_yr_adM <- s_yr_ad 
  
  ## Adults (with sex difference)
  
  # s_yr_saF ~ dunif(0.50, 0.70) 
  # s_yr_adF ~ dunif(0.60, 0.80)  
  # s_yr_saM ~ dunif(0.50, 0.70) 
  # s_yr_adM ~ dunif(0.60, 0.80)  
  
  
  sF_NB[2:3] <- sqrt(s_yr_saF) 
  sF_BN[2:3] <- sqrt(s_yr_saF) 
  
  sF_NB[4] <- sqrt(s_yr_adF) 
  sF_BN[4] <- sqrt(s_yr_adF) 
  
  sM_NB[2:3] <- sqrt(s_yr_saM) 
  sM_BN[2:3] <- sqrt(s_yr_saM) 
  
  sM_NB[4] <- sqrt(s_yr_adM) 
  sM_BN[4] <- sqrt(s_yr_adM) 
  
  
  # Breeding Probability
  
  for (y in 1:ymax) {
    
  Rep[y] ~ dbern(pRep) 
    
  }
  
  pRep ~ dunif(0, 1)
  
  
  # Productivity
  
  for (x in 1:xmax) {
    
    BroodSize[x] ~ dpois(rho[Year_BS[x]])
    
  }
  
  for (t in 1:Tmax){
    
    rho[t] <- mean.rho
    
  }
  
  mean.rho ~ dunif(1, 5)

  
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
  
  for (t in 1:Tmax){
    p[t] <- 0.9
  }
  
  # Population count data (state-space model)
  # Model for the initial population size in Breeding season
  
  for(a in 1:4){
    pinit[a] ~ dunif(0, 200)
    NBreedF[a,1] ~ dpois(pinit[a])  # Female
    NBreedM[a,1] ~ dpois(pinit[a])  # Male
  }
  
  # Book-keeping
  
  NNonF[1:Amax,1] <- 0   # Female
  
  NNonM[1:Amax,1] <- 0   # Male
  
  
  # Process model: Breeding -> Non-Breeding season transition
  
  for (t in 1:Tmax){
    
    # Total number of chicks
    
    Fec[t] ~ dpois(sum(NBreedF[3:4,t]) * pRep * rho[t])
    
    
    # Allocate chicks to a sex
    
    NNonF[1,t+1] ~ dbin(gamma, Fec[t]) # Female chicks 
    NNonM[1,t+1] <- Fec[t] - NNonF[1,t+1] # Male chicks
    
    
    # Survival
    
    for(a in 2:3){
      NNonF[a,t+1] ~ dbin(sF_BN[a-1], NBreedF[a-1,t]) # Female
      NNonM[a,t+1] ~ dbin(sM_BN[a-1], NBreedM[a-1,t]) # Male
    }
    
    # Female
    
    NNonF[4,t+1] <- surv_NBreedF3[t+1] + surv_NBreedF4[t+1] 
    surv_NBreedF3[t+1] ~ dbin(sF_BN[3], NBreedF[3,t])
    surv_NBreedF4[t+1] ~ dbin(sF_BN[4], NBreedF[4,t])
    
    # Male
    
    NNonM[4,t+1] <- surv_NBreedM3[t+1] + surv_NBreedM4[t+1]
    surv_NBreedM3[t+1] ~ dbin(sM_BN[3], NBreedM[3,t])
    surv_NBreedM4[t+1] ~ dbin(sM_BN[4], NBreedM[4,t])
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:4){
      NBreedF[a,t+1] ~ dbin(sF_NB[a], NNonF[a,t+1])  # Female
      NBreedM[a,t+1] ~ dbin(sM_NB[a], NNonM[a,t+1])  # Male
      
    }
  }   
  
  # Observation Model in Non-Breeding
  
  for(j in 1:12){
    
    # Female
    
    ChF_NB[j] ~ dpois(p[NB_yr[j]] * NNonF[1,NB_yr[j]])
    AF_NB[j] ~ dpois(p[NB_yr[j]] * sum(NNonF[2:4,NB_yr[j]]))
    
    # Male 
    
    for(a in 1:4){
      M_NB[a,j] ~ dpois(p[NB_yr[j]] * NNonM[a,NB_yr[j]])
    }
  }
  
  # Observation Model in Breeding
  
  for(h in 1:14){
    
    # Female
    
    JuF_BN[h] ~ dpois(p[BN_yr[h]] * NBreedF[1,BN_yr[h]])
    AF_BN[h] ~ dpois(p[BN_yr[h]] * sum(NBreedF[2:4,BN_yr[h]]))
    
    # Male
    for(a in 1:4){
      M_BN[a,h] ~ dpois(p[BN_yr[h]] * NBreedM[a,BN_yr[h]])
    } 
  }
  
}
)

# Initial values

source("R/GPeafowlIPM_InitialSim_TwoSex_Matrix_BreedProb.R")


Inits <- GP_IPM_Init(Tmax = ny.data + ny.sim, mean.p = 0.9, constant_p = TRUE,
                     survSexDiff = FALSE)
Inits


# Parameters monitored
parameters <- c("sF_NB", "sF_BN","sM_NB", "sM_BN", "mean.rho","gamma",
                "rho", "p", "NBreedF", "NBreedM", "NNonF", "NNonM", 
                "Fec")


# MCMC settings

# ni <- 10
# nb <- 0
# nt <- 1
# nc <- 3

# ni <- 200000     # Run Time around 15 minutes
# nb <- 50000      # Using initial values
# nt <- 30
# nc <- 4

# ni <- 500000     # Run Time around 30 minutes
# nb <- 1000      # Using initial values
# nt <- 30
# nc <- 4

ni <- 10000
nb <- 5000
nt <- 1
nc <- 3



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

# Save output
saveRDS(out, file = "PeafowlIPM_TwoSex_Matrix_TestRun_BreedProb.rds")

plot(out, ask = T)

print(out,2)
summary(out)

# MCMC summaries
print(out,2)

library(MCMCvis)

MCMCsummary(out, params = 'NBreedF', round = 2)
MCMCsummary(out, params = 'NBreedM', round = 2)

MCMCsummary(out, params = 'NNonF', round = 2)
MCMCsummary(out, params = 'NNonM', round = 2)

MCMCsummary(out, params = 'Fec', round = 2)

MCMCsummary(out, params = 'rho', round = 2)

MCMCsummary(out, params = 'p', round = 2)

MCMCsummary(out, 
            params = c('s_NB','s_BN','gamma',
                       'mean.rho', 'sigma.rho'#,
                       #'mean.p','sigma.p'
            ), round = 2)

# Trace Plot

MCMCtrace(out, params = 'NBreed',
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)

MCMCtrace(out, params = 'NNon',
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
