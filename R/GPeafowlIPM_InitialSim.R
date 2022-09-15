
GP_IPM_Init <- function(Tmax){
  Amax <- 4                      # set Tmax and Amax as constants
  
  
  ## Set up vectors and matrices
  NBreed <- NNon <- matrix(NA, nrow = Amax, ncol = Tmax+1)
  s_NB <- s_BN <- rep(NA, Amax)
  pinit <- rep(NA, Amax)
  
  Fec <- rep(NA, Tmax)
  p <- logit.p <- rep(NA, Tmax) 
  rho <- log.rho <- rep(NA, Tmax)
  
  surv_NBreed3 <- surv_NBreed4 <- rep(NA, Tmax+1)
  
  
  ## Sample values for parameters with priors
  # Productivity
  mean.rho <- runif(1, 3, 5)
  sigma.rho <- runif(1, 0, 0.1)
  gamma <- 0.5
  
  
  # Detection Probability
  mean.p <- runif(1, 0.4, 1)
  sigma.p <- runif(1, 0, 1)
  
  
  # Initial Population Sizes
  
  # pinit[1:Amax] <- runif(Amax, 30, 200)
  for(a in 1:Amax){
    pinit[a] <- runif(1, 100, 200)
    NBreed[a,1] <- rpois(1, pinit[a])
  }
  
  NNon[1:Amax,1] <- 0
  surv_NBreed3[1] <- surv_NBreed4[1] <- 0
  
  # Survival rates
  s_NB[1] <- runif(1, 0.30, 0.40)
  s_BN[1] <- runif(1, 0.60, 0.80)
  
  s_yr_sa <- runif(1, 0.50, 0.70) # sub adult (whole year)
  s_NB[2:3] <- sqrt(s_yr_sa) # sub adult (half year, N -> B)
  s_BN[2:3] <- sqrt(s_yr_sa) # sub adult (half year, B -> N)
  
  s_yr_ad <- runif(1, 0.60, 0.80)  # breeder (whole year)
  s_NB[4] <- sqrt(s_yr_ad) # breeder (half year, N -> B)
  s_BN[4] <- sqrt(s_yr_ad) # breeder (half year, B -> N)
  
  
  
  # Calculate vital rates
  
  for (t in 1:Tmax){
    logit.p[t] <- rnorm(1, mean.p, sigma.p)
    p[t] <- plogis(logit.p[t])
  }  
  
  for (t in 1:Tmax){
    log.rho[t] <- rnorm(1, log(mean.rho), sigma.rho)
    rho[t] <- exp(log.rho[t])
  }
  
  #----------------------------------------------------------
  
  # Project population size
  
  for(t in 1:Tmax){
    
    
    # Process model: Breeding -> Non-Breeding season transition    
    
    # Total number of chicks
    Fec[t] <- rpois(1, sum(NBreed[3:4,t]) * rho[t])
    
    # Allocate chicks to a sex
    NNon[1,t+1] <- rbinom(1, Fec[t], gamma) # Female chicks 
    
    # Survival
    for(a in 2:3){
      NNon[a,t+1] <- rbinom(1, NBreed[a-1,t], s_BN[a-1])
    }
    
    surv_NBreed3[t+1] <- rbinom(1, NBreed[3,t], s_BN[3])
    surv_NBreed4[t+1] <- rbinom(1, NBreed[4,t], s_BN[4])
    
    NNon[4,t+1] <- surv_NBreed3[t+1] + surv_NBreed4[t+1]
    
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:Amax){
      NBreed[a,t+1] <- rbinom(1, NNon[a,t+1], s_NB[a])
    }
    
    
  }
  
  
  
  # Arrange as list and return
  
  Inits <- list(
    
    NBreed = NBreed,
    NNon = NNon,
    s_NB = s_NB,
    s_BN = s_BN,
    s_yr_sa = s_yr_sa,
    s_yr_ad = s_yr_ad,
    pinit = pinit,
    Fec = Fec,
    p = p,
    logit.p = logit.p,
    mean.p = mean.p,
    sigma.p = sigma.p,
    rho = rho,
    log.rho = log.rho,
    mean.rho = mean.rho,
    sigma.rho = sigma.rho,
    gamma = gamma,
    surv_NBreed3 = surv_NBreed3,
    surv_NBreed4 = surv_NBreed4
    
  )
  
  Inits
  
  return(Inits)
  
  
}

