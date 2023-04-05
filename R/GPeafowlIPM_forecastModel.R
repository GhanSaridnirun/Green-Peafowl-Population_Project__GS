
library(nimble)

mySeed <- 0
set.seed(mySeed)

# Switch for estimation of rho / use of clutch size data in process model
estimate.rho <- FALSE
# NOTE: The code as it is now only works with estimate.rho = FALSE

## Data Bundle

# Count in Non-Breeding Season

# Year Label
NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2) + 1 

ChF_NB <- c(16,26,69,80,60,20,25,24,48,59,11,12)    #Female Chicks Count

AF_NB <- c(71,69,134,108,86,33,101,68,71,90,34,13)  #All Female 

# Count in Breeding Season

# Year Label

BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)

JuF_BN <- c(9,18,33,33,10,13,6,9,13,27,14,25,11,0)   #Female Juvenile count

AF_BN <- c(18,49,79,71,47,44,30,25,38,54,48,97,17,1) #All Female


#  Male count Data

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
ny.sim <- 20 # Number of years to simulate after the data collection

# Reproduction data

source("R/ReproductionDataPrep.R")

## Set up forecasts

# Define perturbation factor
pert.fac <- 0.1 # Defined as "proportion increase/decrease"

# Set change year (year in which we start a "treatment")
t.change <- ny.data + 8

# List names of vital rates that may be perturbed
VR.names <- c("sF_NB", # Female chick survival
              "sF_BN", # Female juvenile survival
              "sM_NB", # Male chick survival (= female chick survival)
              "sM_BN", # Male juvenile survival (= female juvenile survival)
              "s_yr_saF", # Female subadult annual survival
              "s_yr_saM", # Male subadult annual survival
              "s_yr_adF", # Female adult annual survival
              "s_yr_adM", # Male adult annual survival
              "pRep", # Breeding probability
              "CS", # Clutch size
              "S_C" # Clutch survival
              )

# Set up perturbation matrix
# This matrix has one entry per vital rate per year. 
VR.pert <- matrix(1, nrow = length(VR.names), ncol = ny.data + ny.sim,
                  dimnames = list(VR.names, NULL))

# Apply perturbations as desired
list.VRs_to_perturb <- c("sF_BN", "sM_BN", "s_yr_saF", "s_yr_saM", "s_yr_adF", "s_yr_adM")
# --> In this example, we will apply the perturbation to survival of all juveniles and older birds
for(i in list.VRs_to_perturb){
  VR.pert[which(VR.names == i), t.change:ncol(VR.pert)] <- 1 + pert.fac
}


## Arrange constants
GP.IPMconstants <- list(Tmax = ny.data + ny.sim, 
                        Amax = 4,
                        ny.sim = ny.sim,
                        ny.data = ny.data,
                        NB_yr = NB_yr,
                        BN_yr = BN_yr,
                        xmax = xmax,
                        ymax = ymax,
                        Year_BS = Year_BS,
                        cmax = cmax,
                        VR.pert = VR.pert)

## Arrange data
GP.IPMdata <- list(ChF_NB = ChF_NB,
                   JuF_BN = JuF_BN,
                   AF_NB = AF_NB,
                   AF_BN = AF_BN,
                   M_BN = M_BN,
                   M_NB = M_NB,
                   Rep = Rep,
                   BroodSize = BroodSize,
                   ClutchSize = ClutchSize) 


GPIPM <- list(GP.IPMconstants, GP.IPMdata)  
str(GPIPM)  





## Nimble code for run the whole model

GP.IPMcode <- nimbleCode({
  
  
  # Priors and linear models
  
  # Survival 
  
  ## Chicks and juveniles
  
  Mu.sChick ~ dunif(0.50, 0.60) 
  Mu.sJuv ~ dunif(0.50, 0.60) 
  
  sF_NB[1, 1:Tmax] <- Mu.sChick * VR.pert[1, 1:Tmax]
  sF_BN[1, 1:Tmax] <- Mu.sJuv * VR.pert[2, 1:Tmax] 
  
  sM_NB[1, 1:Tmax] <- Mu.sChick * VR.pert[3, 1:Tmax] 
  sM_BN[1, 1:Tmax] <- Mu.sJuv * VR.pert[4, 1:Tmax]  
  
  ## Adults (no sex difference)
  
  s_yr_sa ~ dunif(0.80, 0.90) 
  s_yr_ad ~ dunif(0.80, 0.90)  
  
  s_yr_saF[1:Tmax] <- s_yr_sa * VR.pert[5,1:Tmax]
  s_yr_saM[1:Tmax] <- s_yr_sa * VR.pert[6,1:Tmax] 
  
  s_yr_adF[1:Tmax] <- s_yr_ad * VR.pert[7,1:Tmax]
  s_yr_adM[1:Tmax] <- s_yr_ad * VR.pert[8,1:Tmax]
  
  for(t in 1:Tmax){
    sF_NB[2:3, t] <- sqrt(s_yr_saF[t]) 
    sF_BN[2:3, t] <- sqrt(s_yr_saF[t]) 
    
    sF_NB[4, t] <- sqrt(s_yr_adF[t]) 
    sF_BN[4, t] <- sqrt(s_yr_adF[t]) 
    
    sM_NB[2:3, t] <- sqrt(s_yr_saM[t]) 
    sM_BN[2:3, t] <- sqrt(s_yr_saM[t]) 
    
    sM_NB[4, t] <- sqrt(s_yr_adM[t]) 
    sM_BN[4, t] <- sqrt(s_yr_adM[t]) 
  }
  
  # Breeding Probability
  
  for (y in 1:ymax) {
    
    Rep[y] ~ dbern(Mu.pRep) 
    
  }
  
  Mu.pRep ~ dunif(0, 1)
  pRep[1:Tmax] <- Mu.pRep * VR.pert[9, 1:Tmax]
  
  # Brood Size and clutch survival [egg to chick]
  
  for (x in 1:xmax) {
    
    BroodSize[x] ~ dpois(rho[Year_BS[x]])
    
  }
  
  if(estimate.rho){
    
    rho[1:Tmax] <- mean.rho
    S_C[1:Tmax] <- (rho[1:Tmax]/mean.CS) 
    
  }else{
    
    rho[1:Tmax] <- CS[1:Tmax] * S_C[1:Tmax]
    S_C[1:Tmax] <- mean.S_C *VR.pert[11, 1:Tmax] 
  }
  
  
  mean.rho ~ dunif(1, 5)
  mean.S_C ~ dunif(0, 1)
  
  
  # Clutch Size
  
  for (x in 1:cmax) {
    
    ClutchSize[x] ~ dpois(mean.CS)
    
  }
  
  mean.CS ~ dunif(3, 11)
  CS[1:Tmax] <- mean.CS * VR.pert[10, 1:Tmax]
  
  
  # Sex ratio of the chicks
  
  gamma <- 0.5
  
  
  
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
    
    Fec[t] ~ dpois(sum(NBreedF[3:4,t]) * pRep[t] * rho[t])
    
    
    # Allocate chicks to a sex
    
    NNonF[1,t+1] ~ dbin(gamma, Fec[t]) # Female chicks 
    NNonM[1,t+1] <- Fec[t] - NNonF[1,t+1] # Male chicks
    
    
    # Survival
    
    for(a in 2:3){
      NNonF[a,t+1] ~ dbin(sF_BN[a-1,t], NBreedF[a-1,t]) # Female
      NNonM[a,t+1] ~ dbin(sM_BN[a-1,t], NBreedM[a-1,t]) # Male
    }
    
    # Female
    
    NNonF[4,t+1] <- surv_NBreedF3[t+1] + surv_NBreedF4[t+1] 
    surv_NBreedF3[t+1] ~ dbin(sF_BN[3,t], NBreedF[3,t])
    surv_NBreedF4[t+1] ~ dbin(sF_BN[4,t], NBreedF[4,t])
    
    # Male
    
    NNonM[4,t+1] <- surv_NBreedM3[t+1] + surv_NBreedM4[t+1]
    surv_NBreedM3[t+1] ~ dbin(sM_BN[3,t], NBreedM[3,t])
    surv_NBreedM4[t+1] ~ dbin(sM_BN[4,t], NBreedM[4,t])
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:4){
      NBreedF[a,t+1] ~ dbin(sF_NB[a,t], NNonF[a,t+1])  # Female
      NBreedM[a,t+1] ~ dbin(sM_NB[a,t], NNonM[a,t+1])  # Male
      
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
# TODO: Make a copy of initial value function and update so that it works with this new model structure
source("R/GPeafowlIPM_InitialSim_TwoSex_Matrix_BreedProb.R")

Inits <- GP_IPM_Init_Pert(Tmax = ny.data + ny.sim, mean.p = 0.9, constant_p = TRUE,
                     survSexDiff = FALSE)
Inits


# Parameters monitored
# TODO: Update parameters to monitor
parameters <- c("sF_NB", "sF_BN","sM_NB", "sM_BN", 
                "NBreedF", "NBreedM", "NNonF", "NNonM", 
                "mean.rho", "rho", "mean.S_C", "S_C", "mean.CS",
                "Fec", "Mu.pRep")


# MCMC settings

# ni <- 10
# nb <- 0
# nt <- 1
# nc <- 3

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
if(estimate.rho){
  saveRDS(out, file = "GPeafowlIPM_forecastModel_rhoEst.rds")
}else{
  saveRDS(out, file = "GPeafowlIPM_forecastModel_rhoDeriv.rds")
}
