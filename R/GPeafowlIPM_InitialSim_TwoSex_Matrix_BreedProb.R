

GP_IPM_Init <- function(Tmax, mean.p, constant_p, survSexDiff){
  
  Amax <- 4                      # set Tmax and Amax as constants
  
  
  ## Set up vectors and matrices
  NBreedF <- NNonF <- matrix(NA, nrow = Amax, ncol = Tmax+1)
  NBreedM <- NNonM <- matrix(NA, nrow = Amax, ncol = Tmax+1)
  sM_NB <- sM_BN <- sF_NB <- sF_BN <- rep(NA, Amax)
  pinit <- rep(NA, Amax)
  Fec <- rep(NA, Tmax)
  p <- logit.p <- rep(NA, Tmax) 
  rho <- rep(NA, Tmax)
  S_C <- rep(NA, Tmax)
  
  surv_NBreedF3 <- surv_NBreedF4 <- rep(NA, Tmax+1)
  surv_NBreedM3 <- surv_NBreedM4 <- rep(NA, Tmax+1)
  
  
  
  ## Sample values for parameters with priors
  
  # Breeding Probability
  
  pRep <- runif(1, 0, 1)
  
  
  # Productivity
  
  
  # Brood Size
  
  mean.rho <- runif(1, 1, 5)
  gamma <- 0.5
  

  # Clutch Size
  
  mean.CS <- runif(1, 3, 11) 
  
  
  # Detection Probability
  
  if(constant_p){
    mean.p <- mean.p
    sigma.p <- 0
  }else{
    mean.p <- runif(1, 0.9*mean.p, 1.1*mean.p)
    sigma.p <- runif(1, 0, 1)
  }
  
  
  # Initial Population Sizes
  
  # pinit[1:Amax] <- runif(Amax, 30, 200)
  for(a in 1:Amax){
    pinit[a] <- runif(1, 0, 200)
    
    NBreedF[a,1] <- rpois(1, pinit[a])  # Female
    NBreedM[a,1] <- rpois(1, pinit[a])  # Male
  }
  
  # Female
  
  NNonF[1:Amax,1] <- 0
  surv_NBreedF3[1] <- surv_NBreedF4[1] <- 0
  
  # Male
  
  NNonM[1:Amax,1] <- 0
  surv_NBreedM3[1] <- surv_NBreedM4[1] <- 0
  
  # Survival rates
  
  ## Chicks and juveniles
  
  sF_NB[1] <- runif(1, 0.50, 0.60) 
  sF_BN[1] <- runif(1, 0.50, 0.60) 
  
  sM_NB[1] <- sF_NB[1] 
  sM_BN[1] <- sF_BN[1]
  
  ## Sub-adult and Adult
  
  if(survSexDiff) {
    
    # ## Adults (with sex difference)------------(TRUE)
    
    s_yr_saF <- runif(1, 0.50, 0.70)
    s_yr_adF <- runif(1, 0.60, 0.80)
    s_yr_saM <- runif(1, 0.50, 0.70)
    s_yr_adM <- runif(1, 0.60, 0.80)
    
  } else {
    
    #   ## Adults (no sex difference)--------------(FALSE)
    
    s_yr_sa <- runif(1, 0.80, 0.90)
    s_yr_ad <- runif(1, 0.80, 0.90)
    
    s_yr_saF <- s_yr_sa
    s_yr_saM <- s_yr_sa
    s_yr_adF <- s_yr_ad
    s_yr_adM <- s_yr_ad
    
  }
  
  sF_NB[2:3] <- sqrt(s_yr_saF)
  sF_BN[2:3] <- sqrt(s_yr_saF)
  sF_NB[4] <- sqrt(s_yr_adF)
  sF_BN[4] <- sqrt(s_yr_adF)
  
  sM_NB[2:3] <- sqrt(s_yr_saM)
  sM_BN[2:3] <- sqrt(s_yr_saM)
  sM_NB[4] <- sqrt(s_yr_adM)
  sM_BN[4] <- sqrt(s_yr_adM)
  
  
  # Calculate vital rates
  
  for (t in 1:Tmax){
    logit.p[t] <- rnorm(1, qlogis(mean.p), sigma.p)
    p[t] <- plogis(logit.p[t])
  }
  
  
  for (t in 1:Tmax){
    rho[t] <- mean.rho
  }
  
  
  # Clutch Survival [from egg to chick]
  
   S_C[1:Tmax] <- rho[1:Tmax]/mean.CS 

   mean.S_C <- mean.rho/mean.CS
   
   
  #----------------------------------------------------------
  
  # Project population size
  
  for(t in 1:Tmax){
    
    
    # Process model: Breeding -> Non-Breeding season transition    
    
    # Total number of chicks
    Fec[t] <- rpois(1, sum(NBreedF[3:4,t]) * pRep * mean.CS * S_C[t])
    
    # Allocate chicks to a sex
    NNonF[1,t+1] <- rbinom(1, Fec[t], gamma) # Female chicks 
    NNonM[1,t+1] <- Fec[t] - NNonF[1,t+1] # Male chicks 
    
    # Survival
    for(a in 2:3){
      NNonF[a,t+1] <- rbinom(1, NBreedF[a-1,t], sF_BN[a-1])  # Female
      NNonM[a,t+1] <- rbinom(1, NBreedM[a-1,t], sM_BN[a-1])  # Male
    }
    
    # Female
    
    surv_NBreedF3[t+1] <- rbinom(1, NBreedF[3,t], sF_BN[3])
    surv_NBreedF4[t+1] <- rbinom(1, NBreedF[4,t], sF_BN[4])
    
    NNonF[4,t+1] <- surv_NBreedF3[t+1] + surv_NBreedF4[t+1]
    
    # Male
    
    surv_NBreedM3[t+1] <- rbinom(1, NBreedM[3,t], sM_BN[3])
    surv_NBreedM4[t+1] <- rbinom(1, NBreedM[4,t], sM_BN[4])
    
    NNonM[4,t+1] <- surv_NBreedM3[t+1] + surv_NBreedM4[t+1]
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:Amax){
      NBreedF[a,t+1] <- rbinom(1, NNonF[a,t+1], sF_NB[a])  #
      NBreedM[a,t+1] <- rbinom(1, NNonM[a,t+1], sM_NB[a])
    }
    
    
  }
  
  
  
  # Arrange as list and return
  
  Inits <- list(  
    NBreedF = NBreedF,
    NBreedM = NBreedM,
    NNonF = NNonF,
    NNonM = NNonM,
    sF_NB = sF_NB,
    sF_BN = sF_BN,
    sM_NB = sM_NB,
    sM_BN = sM_BN,
    s_yr_saF = s_yr_saF,
    s_yr_adF = s_yr_adF,
    s_yr_saM = s_yr_saM,
    s_yr_adM = s_yr_adM,
    mean.CS = mean.CS,
    S_C = S_C,
    mean.S_C = mean.S_C,
    pRep = pRep,
    pinit = pinit,
    Fec = Fec,
    p = p,
    # logit.p = logit.p,
    # mean.p = mean.p,
    # sigma.p = sigma.p,
    rho = rho,
    mean.rho = mean.rho,
    # gamma = gamma,
    surv_NBreedF3 = surv_NBreedF3,
    surv_NBreedF4 = surv_NBreedF4,
    surv_NBreedM3 = surv_NBreedM3,
    surv_NBreedM4 = surv_NBreedM4    
  )
  
  if(!survSexDiff){
    Inits$s_yr_sa <- s_yr_sa
    Inits$s_yr_ad <- s_yr_ad 
  }
  
  Inits
  
  return(Inits)
  
}


#---------------------------------------------------------------


##
##  Update the initial values function for using with GPeafowlIPM_forecastModel.R
##



GP_IPM_Init_Pert <- function(Tmax, VR.pert, mean.p, constant_p, survSexDiff){
  
  Amax <- 4                      # set Tmax and Amax as constants
  
  ## Set up vectors and matrices
  NBreedF <- NNonF <- matrix(NA, nrow = Amax, ncol = Tmax+1)
  NBreedM <- NNonM <- matrix(NA, nrow = Amax, ncol = Tmax+1)
  sM_NB <- sM_BN <- sF_NB <- sF_BN <- matrix(NA, nrow = Amax, ncol = Tmax)
  pinit <- rep(NA, Amax)
  Fec <- rep(NA, Tmax)
  pRep <- rep(NA, Tmax)
  CS <- rep(NA, Tmax)
  p <- logit.p <- rep(NA, Tmax) 
  rho <- rep(NA, Tmax)
  S_C <- rep(NA, Tmax)
  
  surv_NBreedF3 <- surv_NBreedF4 <- rep(NA, Tmax+1)
  surv_NBreedM3 <- surv_NBreedM4 <- rep(NA, Tmax+1)
  
  
  
  ## Sample values for parameters with priors
  
  # Breeding Probability
  
  Mu.pRep <- runif(1, 0, 1)
  
  # Productivity
  # Brood Size
  
  mean.rho <- runif(1, 1, 5)
  gamma <- 0.5
  
  
  # Clutch Size
  
  mean.CS <- runif(1, 3, 11) 
  
  
  # Clutch survival
  
  mean.S_C <- runif(1, 0, 1)
  
  
  # Detection Probability
  
  if(constant_p){
    mean.p <- mean.p
    sigma.p <- 0
  }else{
    mean.p <- runif(1, 0.9*mean.p, 1.1*mean.p)
    sigma.p <- runif(1, 0, 1)
  }
  
  
  # Initial Population Sizes
  
  # pinit[1:Amax] <- runif(Amax, 30, 200)
  for(a in 1:Amax){
    pinit[a] <- runif(1, 0, 200)
    
    NBreedF[a,1] <- rpois(1, pinit[a])  # Female
    NBreedM[a,1] <- rpois(1, pinit[a])  # Male
  }
  
  # Female
  
  NNonF[1:Amax,1] <- 0
  surv_NBreedF3[1] <- surv_NBreedF4[1] <- 0
  
  # Male
  
  NNonM[1:Amax,1] <- 0
  surv_NBreedM3[1] <- surv_NBreedM4[1] <- 0
  
  # Survival rates
  
  ## Chicks and juveniles
  
  Mu.sChick <- runif(1, 0.50, 0.60) 
  Mu.sJuv <- runif(1, 0.50, 0.60) 
  
  sF_NB[1, 1:Tmax] <- pertSurv.nimble(Surv = Mu.sChick, pertFac_t = VR.pert[1, 1:Tmax])
  sF_BN[1, 1:Tmax] <- pertSurv.nimble(Surv = Mu.sJuv, pertFac_t = VR.pert[2, 1:Tmax]) 
  
  sM_NB[1, 1:Tmax] <- pertSurv.nimble(Surv = Mu.sChick, pertFac_t = VR.pert[3, 1:Tmax])
  sM_BN[1, 1:Tmax] <- pertSurv.nimble(Surv = Mu.sJuv, pertFac_t = VR.pert[4, 1:Tmax])  
  
  ## Sub-adult and Adult
  if(survSexDiff) {
    
    # ## Adults (with sex difference)------------(TRUE)
    
    s_yr_saF <- runif(1, 0.50, 0.70)
    s_yr_adF <- runif(1, 0.60, 0.80)
    s_yr_saM <- runif(1, 0.50, 0.70)
    s_yr_adM <- runif(1, 0.60, 0.80)
    
  } else {
    
    #   ## Adults (no sex difference)--------------(FALSE)
    
    s_yr_sa <- runif(1, 0.80, 0.90)
    s_yr_ad <- runif(1, 0.80, 0.90)
    
    s_yr_saF <- s_yr_sa
    s_yr_saM <- s_yr_sa
    s_yr_adF <- s_yr_ad
    s_yr_adM <- s_yr_ad
    
  }
  
  s_yr_saF[1:Tmax] <- pertSurv.nimble(Surv = s_yr_ad, pertFac_t = VR.pert[5,1:Tmax])
  s_yr_saM[1:Tmax] <- pertSurv.nimble(Surv = s_yr_ad, pertFac_t = VR.pert[6,1:Tmax])
  
  s_yr_adF[1:Tmax] <- pertSurv.nimble(Surv = s_yr_ad, pertFac_t = VR.pert[7,1:Tmax])
  s_yr_adM[1:Tmax] <- pertSurv.nimble(Surv = s_yr_ad, pertFac_t = VR.pert[8,1:Tmax])
  
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
  
  pRep[1:Tmax] <- Mu.pRep * VR.pert[9, 1:Tmax]
  
  
  # Clutch Size
  
  CS[1:Tmax] <- mean.CS * VR.pert[10, 1:Tmax]
  
  
  # Clutch Survival [from egg to chick]
  
  # S_C[1:Tmax] <- rho[1:Tmax]/mean.CS 
  # mean.S_C <- mean.rho/mean.CS
  
  S_C[1:Tmax] <- mean.S_C * VR.pert[11, 1:Tmax] 
  rho[1:Tmax] <- CS[1:Tmax] * S_C[1:Tmax]
  
  
  # Observation
  
  for (t in 1:Tmax){
    logit.p[t] <- rnorm(1, qlogis(mean.p), sigma.p)
    p[t] <- plogis(logit.p[t])
  }
  
  
  
  
  #----------------------------------------------------------
  
  # Project population size
  
  for(t in 1:Tmax){
    
    
    # Process model: Breeding -> Non-Breeding season transition    
    
    # Total number of chicks
    Fec[t] <- rpois(1, sum(NBreedF[3:4,t]) * pRep[t] * rho[t])
    
    # Allocate chicks to a sex
    NNonF[1,t+1] <- rbinom(1, Fec[t], gamma) # Female chicks 
    NNonM[1,t+1] <- Fec[t] - NNonF[1,t+1] # Male chicks 
    
    # Survival
    for(a in 2:3){
      NNonF[a,t+1] <- rbinom(1, NBreedF[a-1,t], sF_BN[a-1,t])  # Female
      NNonM[a,t+1] <- rbinom(1, NBreedM[a-1,t], sM_BN[a-1,t])  # Male
    }
    
    # Female
    
    surv_NBreedF3[t+1] <- rbinom(1, NBreedF[3,t], sF_BN[3,t])
    surv_NBreedF4[t+1] <- rbinom(1, NBreedF[4,t], sF_BN[4,t])
    
    NNonF[4,t+1] <- surv_NBreedF3[t+1] + surv_NBreedF4[t+1]
    
    # Male
    
    surv_NBreedM3[t+1] <- rbinom(1, NBreedM[3,t], sM_BN[3,t])
    surv_NBreedM4[t+1] <- rbinom(1, NBreedM[4,t], sM_BN[4,t])
    
    NNonM[4,t+1] <- surv_NBreedM3[t+1] + surv_NBreedM4[t+1]
    
    
    # Process model: Non-Breeding -> Breeding season transition
    
    for(a in 1:Amax){
      NBreedF[a,t+1] <- rbinom(1, NNonF[a,t+1], sF_NB[a,t]) # should be sF_NB[a,t+1])  
      NBreedM[a,t+1] <- rbinom(1, NNonM[a,t+1], sM_NB[a,t]) # should be sM_NB[a,t+1]) 
    }
    
    
  }
  
  
  
  # Arrange as list and return
  
  Inits <- list(  
    NBreedF = NBreedF,
    NBreedM = NBreedM,
    NNonF = NNonF,
    NNonM = NNonM,
    Mu.sChick = Mu.sChick,
    Mu.sJuv = Mu.sJuv,
    sF_NB = sF_NB,
    sF_BN = sF_BN,
    sM_NB = sM_NB,
    sM_BN = sM_BN,
    s_yr_saF = s_yr_saF,
    s_yr_adF = s_yr_adF,
    s_yr_saM = s_yr_saM,
    s_yr_adM = s_yr_adM,
    mean.CS = mean.CS,
    CS = CS,
    S_C = S_C,
    mean.S_C = mean.S_C,
    Mu.pRep = Mu.pRep,
    pRep = pRep,
    pinit = pinit,
    Fec = Fec,
    # p = p,
    # logit.p = logit.p,
    # mean.p = mean.p,
    # sigma.p = sigma.p,
    rho = rho,
    mean.rho = mean.rho,
    # gamma = gamma,
    surv_NBreedF3 = surv_NBreedF3,
    surv_NBreedF4 = surv_NBreedF4,
    surv_NBreedM3 = surv_NBreedM3,
    surv_NBreedM4 = surv_NBreedM4    
  )
  
  if(!survSexDiff){
    Inits$s_yr_sa <- s_yr_sa
    Inits$s_yr_ad <- s_yr_ad 
  }
  
  Inits
  
  return(Inits)
  
}
  
  