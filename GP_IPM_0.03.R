library(nimble)

mySeed <- 0
set.seed(mySeed)


# Data Bundle

#Month order and season for data

month <- factor(c("Nov19","Dec19",
                  "Jan20","Feb20","Mar20","Apr20","May20","Jun20","Jul20","Aug20",
                  "Sep20","Oct20","Nov20","Dec20","Jan21","Feb21","Mar21","Apr21",
                  "May21","Jun21","Jul21","Aug21","Sep21","Oct21","Nov21","Dec21"
),ordered = TRUE)

season <- factor(c("Breeding","Breeding","Breeding","Breeding","Breeding",
                   "Breeding","Non","Non","Non","Non","Non","Non","Breeding","Breeding",
                   "Breeding","Breeding","Breeding","Breeding","Non","Non","Non","Non",
                   "Non","Non","Breeding","Breeding"
),ordered = TRUE)


# Count in Non-Breeding Season

ChF_NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2)             #Season Label for chicks count
ChF_NB <- c(16,26,69,80,60,20,25,24,48,59,11,12)    #Female Chicks Count

Br_NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2)              #Season Label for breeder count
Br_NB <- c(35,34,83,78,56,25,29,28,47,57,14,8)      #Breeder count in Non-Breeding


# Count in Breeding Season

JuF_BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)          #Season Label for juvenile count
JuF_BN <- c(9,18,33,33,10,13,6,9,13,27,14,25,22,0)   #Juvenile count

Br_BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)           #Season Label for breeder count
Br_BN <- c(9,23,32,25,12,10,14,9,11,20,13,27,0,0)    #Breeder count in Breeding


ny=20  # Length = Number of year following Green Peafowl age

str(gp.data <- list(ChF_NB_yr=ChF_NB_yr, ChF_NB=ChF_NB, Br_NB_yr=Br_NB_yr, Br_NB=Br_NB, 
                    JuF_BN_yr=JuF_BN_yr, JuF_BN=JuF_BN, Br_BN_yr=Br_BN_yr, Br_BN=Br_BN,
                    ny=ny)
)



GP.code <- nimbleCode({
  
  
  # Priors and linear models
  # Survival
  
  s_NB[1] ~ dunif(0.30, 0.40)  # chick (half year)
  s_BN[1] ~ dunif(0.40, 0.60)  # juvenile (half year)
  
  s_yr_sa ~ dunif(0.50, 0.70) # sub adult (whole year)
  s_NB[2:3] <- sqrt(s_yr_sa) # sub adult (half year, N -> B)
  s_BN[2:3] <- sqrt(s_yr_sa) # sub adult (half year, B -> N)
  
  s_yr_ad ~ dunif(0.60, 0.80)  # breeder (whole year)
  s_NB[4] <- sqrt(s_yr_ad) # breeder (half year, N -> B)
  s_BN[4] <- sqrt(s_yr_ad) # breeder (half year, B -> N)
  
  
  # Productivity
  for (t in 1:ny){
    log.rho[t] ~ dnorm(log.mean.rho, tau.rho)
    rho[t] <- exp(log.rho[t])
  }
  
  mean.rho ~ dunif(0, 5)
  log.mean.rho <- log(mean.rho)
  sigma.rho ~ dunif(0, 1)
  tau.rho <- pow(sigma.rho, -2)  
  
  # Sex ratio of the chicks
  gamma ~ dbeta(1, 1)  
  
  
  # Detection Probability
  for (t in 1:ny){
    logit.p[t] ~ dnorm(lmean.p, tau.p)
    p[t] <- ilogit(logit.p[t])
  }  
  
  mean.p ~ dbeta(1, 1)
  lmean.p <- logit(mean.p)
  sigma.p ~ dunif(0, 5)
  tau.p <- pow(sigma.p, -2)
  
  
  # Population count data (state-space model)
  # Model for the initial population size in Breeding season
  
  for(a in 1:4){
    pinit[a] ~ dunif(30, 200)
    NBreed[a,1] ~ dpois(pinit[a]) 
  }
  
  
  # Process model: Breeding -> Non-Breeding season transition
  for (t in 1:ny){
    
    # Total number of chicks
    Fec[t] ~ dpois(sum(NBreed[3:4,t]) * rho[t])
    
    # Allocate chicks to a sex
    NNon[1,t+1] ~ dbin(gamma, Fec[t]) # Female chicks 
    
    # Survival
    for(a in 2:3){
      NNon[a,t+1] ~ dbin(s_BN[a-1], NBreed[a-1,t])
    }
    
    Nnon[4,t+1] <- surv_NBreed3[t+1] + surv_NBreed4[t+1]
    surv_NBreed3[t+1] ~ dbin(s_BN[3], NBreed[3,t])
    surv_NBreed4[t+1] ~ dbin(s_BN[4], NBreed[4,t])
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:4){
      NBreed[a,t+1] ~ dbin(s_NB[a], NNon[a,t+1])
    }
    
    
    # This section I think I just forget to add p[t] in the right side 
    # due to it would provide the estimated probability of detection
    
    # Observation Model in Non-Breeding
    for(j in 1:12) {
      
      ChF_NB[j] ~ dpois(NNon[1,ChF_NB_yr[j]])
      Br_NB[j] ~ dpois(NNon[4,Br_NB_yr[j]])
      
    }
    
    # Observation Model in Breeding
    for(h in 1:14) {
      
      JuF_BN[h] ~ dpois(NBreed[2,JuF_BN_yr[h]])
      Br_BN[h] ~ dpois(NBreed[4,Br_BN_yr[h]])
      
    }
  } 
  
}
)

# Initial values
inits <- function() {
  list(
    NBreed ~ runif(1, 30, 200),
    NBreed[1] ~ runif(1, 30, 200),
    NBreed[2] ~ runif(1, 30, 200),
    NBreed[3] ~ runif(1, 30, 200),
    NBreed[4] ~ runif(1, 30, 200),
    NNon ~ runif(1, 30, 200),
    NNon[1] ~ runif(1, 30, 200),
    NNon[2] ~ runif(1, 30, 200),
    NNon[3] ~ runif(1, 30, 200),
    NNon[4] ~ runif(1, 30, 200),    
    s_NB[1] ~ runif(1, 0.30, 0.40),
    s_BN[1] ~ runif(1, 0.40, 0.60),
    s_NB[2] ~ runif(1, 0.50, 0.70),
    s_BN[2] ~ runif(1, 0.50, 0.70),
    s_NB[3] ~ runif(1, 0.50, 0.70),
    s_BN[3] ~ runif(1, 0.50, 0.70),
    s_NB[4] ~ runif(1, 0.60, 0.80),
    s_BN[4] ~ runif(1, 0.60, 0.80),  
    #rho  = runif(1, 0, 10),
    mean.rho = runif(1, 0, 5),
    Fec = runif(1, 2, 4),
    #log.rho = runif
    sigma.rho = runif(1, 0, 1),
    gamma = runif(1, 0.4, 0.6),
    #p = runif(1, 0.4, 1),
    mean.p = runif(1, 0.4, 1),
    sigma.p = runif(1, 0, 1)
  )  
}  



Inits <- list(inits())


# Parameters monitored
parameters <- c("s_NB", "s_BN", "mean.rho", "gamma", "mean.p", "sigma.rho", "sigma.p", "rho", "p",
                "NBreed", "NNon", "Fec")


# MCMC settings
ni <- 10
nb <- 0
nt <- 1
#nc <- 3
nc <- 1

# Run
out <- nimbleMCMC(code = GP.code,
                  constants = gp.data,
                  inits = Inits,
                  monitors = parameters,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = nb,
                  niter = ni,
                  thin = nt,
                  nchains = nc)

print(out)

