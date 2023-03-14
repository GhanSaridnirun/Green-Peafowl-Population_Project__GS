########################
# GREEN PEAFOWL IPM    #
# SENSITIVITY ANALYSIS #
########################


# 1. Set parameter values #
#-------------------------#

## Productivity 
pRep <- 0.9 # Female breeding probability
mean.CS <- 8 # Average clutch size
gamma <- 0.5 # Sex ratio at birth
S_C <- 0.7 # Clutch survival

## Survival
sF_NB <- c(0.55, 0.65, 0.75, 0.9) # Female survival Non-breeding -> Breeding (chicks, 1yr, 2yr, 3yr)
sF_BN <- c(0.55, 0.65, 0.75, 0.9) # Female survival Breeding -> Non-breeding (juveniles, 1yr, 2yr, 3yr)
sM_NB <- c(0.55, 0.65, 0.75, 0.9) # Male survival Non-breeding -> Breeding (chicks, 1yr, 2yr, 3yr)
sM_BN <- c(0.55, 0.65, 0.75, 0.9) # Male survival Breeding -> Non-breeding (juveniles, 1yr, 2yr, 3yr)


# 2. Build projection matrix #
#----------------------------#

## Function to make projection matrix
source("R/make.GPprojMatrix.R")

## Build annual matrix
mat.ann <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                             sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                             seasonal = FALSE)

## Build seasonal matrix
mat.season <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                seasonal = TRUE)


# 3. Calculate asymptotic population growth rate #
#------------------------------------------------#

## Set initial population to use in simulations
N1.ann <- c(100, 100, # Chicks
            25, 22, # 1 yr olds
            20, 15, # 2yr olds
            50, 40) # 3+yr olds

N1.season <- c(N1.ann[1:2],
               0, 0,
               N1.ann[3:4],
               0, 0,
               N1.ann[5:6],
               0, 0,
               N1.ann[7:8],
               0, 0)
# NOTE: The 0's in the seasonal population are for the stages not present in the first season.
# Since we'll consider an annual census at the start of the breeding season, there will be 0
# individuals of all the non-breeding season stages


## Annual matrix
# Using eigen decomposition
lambda.ann <- eigen(mat.ann$A)$values[1]
lambda.ann
# --> Just above one, meaning population growing by about 1.3% per year (under stable conditions)

# Using simulation
tSim <- 200
N.sim.ann <- matrix(NA, nrow = length(N1.ann), ncol = tSim)
N.sim.ann[, 1] <- N1.ann

for(t in 1:(tSim-1)){
  N.sim.ann[, t+1] <- mat.ann$A %*% N.sim.ann[, t]
}
matplot(t(N.sim.ann), type = "l") # Population by stage
plot(colSums(N.sim.ann), type = "l") # Entire population
plot(log(colSums(N.sim.ann)), type = "l") # Entire population (log)

lambda.ann.sim <- sum(N.sim.ann[,tSim]) / sum(N.sim.ann[,tSim-1])

# Comparison
lambda.ann
lambda.ann.sim
# --> They are basically identical (as they should be)


## Seasonal matrices
# Using eigen decomposition
lambda.BtoN <- eigen(mat.season$A_BtoN)$values[1]
lambda.BtoN 
lambda.NtoB <- eigen(mat.season$A_NtoB)$values[1]
lambda.NtoB
# --> Both are 0 because life cycle is "interrupted" when considering seasons separately

# Using simulation
N.sim.season <- matrix(NA, nrow = length(N1.season), ncol = tSim*2) # Twice as many time steps since we have two seasons
N.sim.season[, 1] <- N1.season

for(t in 1:((tSim*2)-1)){
  
  # Non-breeding to breeding season
  if(t %% 2 == 1){ # If t is odd ...
    N.sim.season[, t+1] <- mat.season$A_NtoB %*% N.sim.season[, t]
  } 
  
  # Breeding to non-breeding season
  if(t %% 2 == 0){ # If t is even ...
    N.sim.season[, t+1] <- mat.season$A_BtoN %*% N.sim.season[, t]
  } 
  
}

matplot(t(N.sim.season), type = "l") # Population by stage (both seasons)
plot(colSums(N.sim.season), type = "l") # Entire population (both seasons)
plot(log(colSums(N.sim.season)), type = "l") # Entire population (log, both seasons)

N.sim.NB <- N.sim.season[, which(1:(tSim*2) %% 2 == 1)]
matplot(t(N.sim.NB), type = "l") # Population by stage (both seasons)
plot(colSums(N.sim.NB), type = "l") # Entire population (both seasons)
plot(log(colSums(N.sim.NB)), type = "l") # Entire population (log, both seasons)

lambda.season.sim <- sum(N.sim.NB[,tSim]) / sum(N.sim.NB[,tSim-1])


## Compare (if matrices are specified correctly, they should have the same lambda)
lambda.ann.sim == lambda.season.sim
plot(colSums(N.sim.ann), type = "l")
lines(colSums(N.sim.ann), col = "red")
# --> They are the same, so our matrices are correct.
# --> Therefore we can safely use the simpler, annual matrix in the following



# 4. Calculate matrix element sensitivities #
#-------------------------------------------#

## Using the popbio package
popbio::eigen.analysis(mat.ann$A)
# --> This approach gives us NA for sensitivities
# --> The reason for that is that our projection matrix is singular, and its inverse (necessary for sensitivity calculations) not defined

## Using perturbation analysis
# NOTE: Sensitivities (and elasticities) can also be calculated using perturbation analysis instead of the analytical approach implemented in popbio
#       The perturbation approach will also work for singular matrices, and is very straightforward to apply also for
#       calculating sensitivities/elasticities with respect to vital rates (below).
#       Additionally, you can use it not just for sensitivities/elasticities of asymptotic population growth rate, but for any metric of your choice.

# Set perturbation factor
dy <- 1e-5

# Set up original and perturbed matrix
A_orig <- mat.ann$A
A_pert <- A_orig

# Perturb target element in matrix
A_pert_SChF <- A_pert
A_pert_SChF[3, 1] <- A_pert_SChF[3, 1] + dy

# Calculate population growth rate for both matrices
lam_orig <- as.numeric(eigen(A_orig)$values[1])
lam_pert_SChF <- as.numeric(eigen(A_pert_SChF)$values[1])

# Calculate sensitivity of population growth rate to target element
sens <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens[3, 1] <- (lam_pert_SChF - lam_orig) / dy 


# TODO: 
# - Generalize above code to calculate sensitivities for all (non-0) matrix elements

#Sensitivity of ChM
A_pert_SChM <- A_pert
A_pert_SChM[4, 2] <- A_pert_SChM[4, 2] + dy
lam_pert_SChM <- as.numeric(eigen(A_pert_SChM)$values[1])
sens[4, 2] <- (lam_pert_SChM - lam_orig) / dy

#Sensitivity of 1yF
A_pert_S1yF <- A_pert
A_pert_S1yF[5, 3] <- A_pert_S1yF[5, 3] + dy
lam_pert_S1yF <- as.numeric(eigen(A_pert_S1yF)$values[1])
sens[5, 3] <- (lam_pert_S1yF - lam_orig) / dy

#Sensitivity of 1yM
A_pert_S1yM <- A_pert
A_pert_S1yM[6, 4] <- A_pert_S1yM[6, 4] + dy
lam_pert_S1yM <- as.numeric(eigen(A_pert_S1yM)$values[1])
sens[6, 4] <- (lam_pert_S1yM - lam_orig) / dy

#Sensitivity of 2yF produce ChF
A_pert_S2yFF <- A_pert
A_pert_S2yFF[1, 5] <- A_pert_S2yFF[1, 5] + dy
lam_pert_S2yFF <- as.numeric(eigen(A_pert_S2yFF)$values[1])
sens[1, 5] <- (lam_pert_S2yFF - lam_orig) / dy

#Sensitivity of 2yF produce ChM
A_pert_S2yFM <- A_pert
A_pert_S2yFM[2, 5] <- A_pert_S2yFM[2, 5] + dy 
lam_pert_S2yFM <- as.numeric(eigen(A_pert_S2yFM)$values[1])
sens[2, 5] <- (lam_pert_S2yFM - lam_orig) / dy

#Sensitivity of 2yF
A_pert_S2yF <- A_pert
A_pert_S2yF[7, 5] <- A_pert_S2yF[7, 5] + dy 
lam_pert_S2yF <- as.numeric(eigen(A_pert_S2yF)$values[1])
sens[7, 5] <- (lam_pert_S2yF - lam_orig) / dy 

#Sensitivity of 2yM
A_pert_S2yM <- A_pert
A_pert_S2yM[8, 6] <- A_pert_S2yM[8, 6] + dy
lam_pert_S2yM <- as.numeric(eigen(A_pert_S2yM)$values[1])
sens[8, 6] <- (lam_pert_S2yM - lam_orig) / dy 

#Sensitivity of 3yF produce ChF
A_pert_S3yFF <- A_pert
A_pert_S3yFF[1, 7] <- A_pert_S3yFF[1, 7] + dy
lam_pert_S3yFF <- as.numeric(eigen(A_pert_S3yFF)$values[1])
sens[1, 7] <- (lam_pert_S3yFF - lam_orig) / dy 

#Sensitivity of 3yF produce ChM
A_pert_S3yFM <- A_pert
A_pert_S3yFM[2, 7] <- A_pert_S3yFM[2, 7] + dy
lam_pert_S3yFM <- as.numeric(eigen(A_pert_S3yFM)$values[1])
sens[2, 7] <- (lam_pert_S3yFM - lam_orig) / dy 

#Sensitivity of 3yF
A_pert_S3yF <- A_pert
A_pert_S3yF[7, 7] <- A_pert_S3yF[7, 7] + dy
lam_pert_S3yF <- as.numeric(eigen(A_pert_S3yF)$values[1])
sens[7, 7] <- (lam_pert_S3yF - lam_orig) / dy 

#Sensitivity of 3yM
A_pert_S3yM <- A_pert
A_pert_S3yM[8, 8] <- A_pert_S3yM[8, 8] + dy
lam_pert_S3yM <- as.numeric(eigen(A_pert_S3yM)$values[1])
sens[8, 8] <- (lam_pert_S3yM - lam_orig) / dy 


sens

row.names(sens) <- (c("ChF", "ChM","1yF", "1yM","2yF", "2yM","3yF", "3yM"))
colnames(sens) <- (c("ChF", "ChM","1yF", "1yM","2yF", "2yM","3yF", "3yM"))

sens



# - Visualize sensitivity matrix (tip: start with fields::image.plot)

fields::image.plot(sens)



# 5. Calculate matrix element elasticities #
#------------------------------------------#

## Using perturbation analysis

# Set perturbation factor
dy <- 1e-5

# Set up original and perturbed matrix
A_orig <- mat.ann$A
A_pert <- A_orig

# Perturb target element in matrix
A_pert_EChF <- A_pert
A_pert_EChF[3, 1] <- A_pert_EChF[3, 1] * (1 + dy) 


# Calculate population growth rate for both matrices
lam_orig <- as.numeric(eigen(A_orig)$values[1])
lam_pert_EChF <- as.numeric(eigen(A_pert_EChF)$values[1])

# Calculate sensitivity of population growth rate to target element
elas <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas[3, 1] <- (lam_pert_EChF - lam_orig) / (lam_orig * dy) 


# TODO: 
# - Generalize above code to calculate sensitivities for all (non-0) matrix elements

# Elasticity of ChM
A_pert_EChM <- A_pert
A_pert_EChM[4, 2] <- A_pert_EChM[4, 2] * (1 + dy)  
lam_pert_EChM <- as.numeric(eigen(A_pert_EChM)$values[1])
elas[4, 2] <- (lam_pert_EChM - lam_orig) / (lam_orig * dy) 

# Elasticity of 1yF
A_pert_E1yF <- A_pert
A_pert_E1yF[5, 3] <- A_pert_E1yF[5, 3] * (1 + dy) 
lam_pert_E1yF <- as.numeric(eigen(A_pert_E1yF)$values[1])
elas[5, 3] <- (lam_pert_E1yF - lam_orig) / (lam_orig * dy) 

# Elasticity of 1yM
A_pert_E1yM <- A_pert
A_pert_E1yM[6, 4] <- A_pert_E1yM[6, 4] * (1 + dy)  
lam_pert_E1yM <- as.numeric(eigen(A_pert_E1yM)$values[1])
elas[6, 4] <- (lam_pert_E1yM - lam_orig) / (lam_orig * dy) 

# Elasticity of 2yF produce ChF
A_pert_E2yFF <- A_pert
A_pert_E2yFF[1, 5] <- A_pert_E2yFF[1, 5] * (1 + dy) 
lam_pert_E2yFF <- as.numeric(eigen(A_pert_E2yFF)$values[1])
elas[1, 5] <- (lam_pert_E2yFF - lam_orig) / (lam_orig * dy) 

# Elasticity of 2yF produce ChM
A_pert_E2yFM <- A_pert
A_pert_E2yFM[2, 5] <- A_pert_E2yFM[2, 5] * (1 + dy)  
lam_pert_E2yFM <- as.numeric(eigen(A_pert_E2yFM)$values[1])
elas[2, 5] <- (lam_pert_E2yFM - lam_orig) / (lam_orig * dy) 

# Elasticity of 2yF
A_pert_E2yF <- A_pert
A_pert_E2yF[7, 5] <- A_pert_E2yF[7, 5] * (1 + dy)   
lam_pert_E2yF <- as.numeric(eigen(A_pert_E2yF)$values[1])
elas[7, 5] <- (lam_pert_E2yF - lam_orig) / (lam_orig * dy) 

# Elasticity of 2yM
A_pert_E2yM <- A_pert
A_pert_E2yM[8, 6] <- A_pert_E2yM[8, 6] * (1 + dy)  
lam_pert_E2yM <- as.numeric(eigen(A_pert_E2yM)$values[1])
elas[8, 6] <- (lam_pert_E2yM - lam_orig) / (lam_orig * dy) 

# Elasticity of 3yF produce ChF
A_pert_E3yFF <- A_pert
A_pert_E3yFF[1, 7] <- A_pert_E3yFF[1, 7] * (1 + dy) 
lam_pert_E3yFF <- as.numeric(eigen(A_pert_E3yFF)$values[1])
elas[1, 7] <- (lam_pert_E3yFF - lam_orig) / (lam_orig * dy) 

# Elasticity of 3yF produce ChM
A_pert_E3yFM <- A_pert
A_pert_E3yFM[2, 7] <- A_pert_E3yFM[2, 7] * (1 + dy) 
lam_pert_E3yFM <- as.numeric(eigen(A_pert_E3yFM)$values[1])
elas[2, 7] <- (lam_pert_E3yFM - lam_orig) / (lam_orig * dy) 

# Elasticity of 3yF
A_pert_E3yF <- A_pert
A_pert_E3yF[7, 7] <- A_pert_E3yF[7, 7] * (1 + dy)  
lam_pert_E3yF <- as.numeric(eigen(A_pert_E3yF)$values[1])
elas[7, 7] <- (lam_pert_E3yF - lam_orig) / (lam_orig * dy) 

# Elasticity of 3yM
A_pert_E3yM <- A_pert
A_pert_E3yM[8, 8] <- A_pert_E3yM[8, 8] * (1 + dy) 
lam_pert_E3yM <- as.numeric(eigen(A_pert_E3yM)$values[1])
elas[8, 8] <- (lam_pert_E3yM - lam_orig) / (lam_orig * dy) 


elas

row.names(elas) <- (c("ChF", "ChM","1yF", "1yM","2yF", "2yM","3yF", "3yM"))
colnames(elas) <- (c("ChF", "ChM","1yF", "1yM","2yF", "2yM","3yF", "3yM"))

elas


# - Visualize elasticity matrix (tip: start with fields::image.plot)

fields::image.plot(elas)



# 6. Calculate vital rate sensitivities #
#---------------------------------------#

# TODO: 
# - Adapt the above code to calculate sensitivities with respect to vital rates (instead of matrix elements)
#   NOTE: Keep in mind that unlike matrix elements, vital rates may appear several times in the matrix!
# - Generalize your code to calculate sensitivities for all vital rates

## Using perturbation analysis

# Set perturbation factor
dy <- 1e-5

# Setting the perturbation for sensitivity of target vital rates
S.pRep <- pRep + dy # pertubate Breeding probability
S.mean.CS <- mean.CS + dy # pertubate Mean clutch size
S.S_C <- S_C + dy # pertubate Clutch survival
S.sF_NB_Ch <- sF_NB + c(dy, 0, 0, 0) # pertubate Female chick Non-breeding survival
S.sF_NB_1yr <- sF_NB + c(0, dy, 0, 0) # pertubate Female 1yr Non-breeding survival
S.sF_NB_2yr <- sF_NB + c(0, 0, dy, 0) # pertubate Female 2yr Non-breeding survival
S.sF_NB_3yr <- sF_NB + c(0, 0, 0, dy) # pertubate Female 3yr Non-breeding survival
S.sF_BN_Ju <- sF_BN + c(dy, 0, 0, 0) # pertubate Female Juvenile Breeding survival
S.sF_BN_1yr <- sF_BN + c(0, dy, 0, 0) # pertubate Female 1yr Breeding survival
S.sF_BN_2yr <- sF_BN + c(0, 0, dy, 0) # pertubate Female 2yr Breeding survival
S.sF_BN_3yr <- sF_BN + c(0, 0, 0, dy) # pertubate Female 3yr Breeding survival
S.sM_NB_Ch <- sM_NB + c(dy, 0, 0, 0) # pertubate Male chick Non-breeding survival
S.sM_NB_1yr <- sM_NB + c(0, dy, 0, 0) # pertubate Male 1yr Non-breeding survival
S.sM_NB_2yr <- sM_NB + c(0, 0, dy, 0) # pertubate Male 2yr Non-breeding survival
S.sM_NB_3yr <- sM_NB + c(0, 0, 0, dy) # pertubate Male 3yr Non-breeding survival
S.sM_BN_Ju <- sM_BN + c(dy, 0, 0, 0) # pertubate Male juvenile Breeding survival
S.sM_BN_1yr <- sM_BN + c(0, dy, 0, 0) # pertubate Male 1yr Breeding survival
S.sM_BN_2yr <- sM_BN + c(0, 0, dy, 0) # pertubate Male 2yr Breeding survival
S.sM_BN_3yr <- sM_BN + c(0, 0, 0, dy) # pertubate Male 3yr Breeding survival

# Origin Matrix
Sen_orig <- mat.ann$A

# Build matrix with pertubated target vital rate for sensitivity

# Breeding probability
mat.Sen.pRep <- make.GPprojMatrix(pRep = S.pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                  sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                  seasonal = FALSE)

# Mean clutch size
mat.Sen.mean.CS <- make.GPprojMatrix(pRep = pRep, mean.CS = S.mean.CS, S_C = S_C, gamma = gamma, 
                                     sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                     seasonal = FALSE)

# Clutch survival
mat.Sen.S_C <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S.S_C, gamma = gamma, 
                                 sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                 seasonal = FALSE)

# Female chick Non-breeding survival
mat.Sen.sF_NB_Ch <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = S.sF_NB_Ch, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                      seasonal = FALSE)

# Female 1yr Non-breeding survival
mat.Sen.sF_NB_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = S.sF_NB_1yr, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 2yr Non-breeding survival
mat.Sen.sF_NB_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = S.sF_NB_2yr, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 3yr Non-breeding survival
mat.Sen.sF_NB_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = S.sF_NB_3yr, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female Juvenile Breeding survival
mat.Sen.sF_BN_Ju <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = sF_NB, sF_BN = S.sF_BN_Ju, sM_NB = sM_NB, sM_BN = sM_BN,
                                      seasonal = FALSE)

# Female 1yr Breeding survival
mat.Sen.sF_BN_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = S.sF_BN_1yr, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 2yr Breeding survival
mat.Sen.sF_BN_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = S.sF_BN_2yr, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 3yr Breeding survival
mat.Sen.sF_BN_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = S.sF_BN_3yr, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male chick Non-breeding survival
mat.Sen.sM_NB_Ch <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = S.sM_NB_Ch, sM_BN = sM_BN,
                                      seasonal = FALSE)

# Male 1yr Non-breeding survival
mat.Sen.sM_NB_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = S.sM_NB_1yr, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male 2yr Non-breeding survival
mat.Sen.sM_NB_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = S.sM_NB_2yr, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male 3yr Non-breeding survival
mat.Sen.sM_NB_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = S.sM_NB_3yr, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male juvenile Breeding survival
mat.Sen.sM_BN_Ju <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = S.sM_BN_Ju,
                                      seasonal = FALSE)

# Male 1yr Breeding survival
mat.Sen.sM_BN_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = S.sM_BN_1yr,
                                       seasonal = FALSE)

# Male 2yr Breeding survival
mat.Sen.sM_BN_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = S.sM_BN_2yr,
                                       seasonal = FALSE)

# Male 3yr Breeding survival
mat.Sen.sM_BN_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = S.sM_BN_3yr,
                                       seasonal = FALSE)


# Select only numerical matrix for calculation
Sen.pRep <- mat.Sen.pRep$A # Breeding probability
Sen.mean.CS <- mat.Sen.mean.CS$A
Sen.S_C <- mat.Sen.S_C$A
Sen.sF_NB_Ch <- mat.Sen.sF_NB_Ch$A
Sen.sF_NB_1yr <- mat.Sen.sF_NB_1yr$A 
Sen.sF_NB_2yr <- mat.Sen.sF_NB_2yr$A
Sen.sF_NB_3yr <- mat.Sen.sF_NB_3yr$A
Sen.sF_BN_Ju <- mat.Sen.sF_BN_Ju$A
Sen.sF_BN_1yr <- mat.Sen.sF_BN_1yr$A
Sen.sF_BN_2yr <- mat.Sen.sF_BN_2yr$A
Sen.sF_BN_3yr <- mat.Sen.sF_BN_3yr$A
Sen.sM_NB_Ch <- mat.Sen.sM_NB_Ch$A
Sen.sM_NB_1yr <- mat.Sen.sM_NB_1yr$A
Sen.sM_NB_2yr <- mat.Sen.sM_NB_2yr$A
Sen.sM_NB_3yr <- mat.Sen.sM_NB_3yr$A
Sen.sM_BN_Ju <- mat.Sen.sM_BN_Ju$A
Sen.sM_BN_1yr <- mat.Sen.sM_BN_1yr$A
Sen.sM_BN_2yr <- mat.Sen.sM_BN_2yr$A
Sen.sM_BN_3yr <- mat.Sen.sM_BN_3yr$A

# Calculate population growth rate for Origin matrix
lam_Sen_orig <- as.numeric(eigen(Sen_orig)$values[1])

# Calculate population growth rate for Origin matrix
lam_Sen.pRep <- as.numeric(eigen(Sen.pRep)$values[1])
lam_Sen.mean.CS <- as.numeric(eigen(Sen.mean.CS)$values[1])
lam_Sen.S_C <- as.numeric(eigen(Sen.S_C)$values[1])
lam_Sen.sF_NB_Ch <- as.numeric(eigen(Sen.sF_NB_Ch)$values[1])
lam_Sen.sF_NB_1yr <- as.numeric(eigen(Sen.sF_NB_1yr)$values[1])
lam_Sen.sF_NB_2yr <- as.numeric(eigen(Sen.sF_NB_2yr)$values[1])
lam_Sen.sF_NB_3yr <- as.numeric(eigen(Sen.sF_NB_3yr)$values[1])
lam_Sen.sF_BN_Ju <- as.numeric(eigen(Sen.sF_BN_Ju)$values[1])
lam_Sen.sF_BN_1yr <- as.numeric(eigen(Sen.sF_BN_1yr)$values[1])
lam_Sen.sF_BN_2yr <- as.numeric(eigen(Sen.sF_BN_2yr)$values[1])
lam_Sen.sF_BN_3yr <- as.numeric(eigen(Sen.sF_BN_3yr)$values[1])
lam_Sen.sM_NB_Ch <- as.numeric(eigen(Sen.sM_NB_Ch)$values[1])
lam_Sen.sM_NB_1yr <- as.numeric(eigen(Sen.sM_NB_1yr)$values[1])
lam_Sen.sM_NB_2yr <- as.numeric(eigen(Sen.sM_NB_2yr)$values[1])
lam_Sen.sM_NB_3yr <- as.numeric(eigen(Sen.sM_NB_3yr)$values[1])
lam_Sen.sM_BN_Ju <- as.numeric(eigen(Sen.sM_BN_Ju)$values[1])
lam_Sen.sM_BN_1yr <- as.numeric(eigen(Sen.sM_BN_1yr)$values[1])
lam_Sen.sM_BN_2yr <- as.numeric(eigen(Sen.sM_BN_2yr)$values[1])
lam_Sen.sM_BN_3yr <- as.numeric(eigen(Sen.sM_BN_3yr)$values[1])

# Calculate sensitivity of population growth rate to target element
Sen_pRep <- (lam_Sen.pRep - lam_Sen_orig) / dy
Sen_mean.CS  <- (lam_Sen.mean.CS - lam_Sen_orig) / dy
Sen_S_C <- (lam_Sen.S_C - lam_Sen_orig) / dy
Sen_sF_NB_Ch <- (lam_Sen.sF_NB_Ch - lam_Sen_orig) / dy
Sen_sF_NB_1yr <- (lam_Sen.sF_NB_1yr - lam_Sen_orig) / dy
Sen_sF_NB_2yr <- (lam_Sen.sF_NB_2yr - lam_Sen_orig) / dy
Sen_sF_NB_3yr <- (lam_Sen.sF_NB_3yr - lam_Sen_orig) / dy
Sen_sF_BN_Ju <- (lam_Sen.sF_BN_Ju - lam_Sen_orig) / dy
Sen_sF_BN_1yr <- (lam_Sen.sF_BN_1yr - lam_Sen_orig) / dy
Sen_sF_BN_2yr <- (lam_Sen.sF_BN_2yr - lam_Sen_orig) / dy
Sen_sF_BN_3yr <- (lam_Sen.sF_BN_3yr - lam_Sen_orig) / dy
Sen_sM_NB_Ch <- (lam_Sen.sM_NB_Ch - lam_Sen_orig) / dy
Sen_sM_NB_1yr <- (lam_Sen.sM_NB_1yr - lam_Sen_orig) / dy
Sen_sM_NB_2yr <- (lam_Sen.sM_NB_2yr - lam_Sen_orig) / dy
Sen_sM_NB_3yr <- (lam_Sen.sM_NB_3yr - lam_Sen_orig) / dy
Sen_sM_BN_Ju <- (lam_Sen.sM_BN_Ju - lam_Sen_orig) / dy
Sen_sM_BN_1yr <- (lam_Sen.sM_BN_1yr - lam_Sen_orig) / dy
Sen_sM_BN_2yr <- (lam_Sen.sM_BN_2yr - lam_Sen_orig) / dy
Sen_sM_BN_3yr <- (lam_Sen.sM_BN_3yr - lam_Sen_orig) / dy

#Gather the results
Sensitivities <- list(c(Sen_pRep, Sen_mean.CS, Sen_S_C, 
                        Sen_sF_NB_Ch, Sen_sF_NB_1yr, Sen_sF_NB_2yr, Sen_sF_NB_3yr,
                        Sen_sF_BN_Ju, Sen_sF_BN_1yr, Sen_sF_BN_2yr, Sen_sF_BN_3yr,
                        Sen_sM_NB_Ch, Sen_sM_NB_1yr, Sen_sM_NB_2yr, Sen_sM_NB_3yr,
                        Sen_sM_BN_Ju, Sen_sM_BN_1yr, Sen_sM_BN_2yr, Sen_sM_BN_3yr))



# - Visualize sensitivities for all vital rates 



# 7. Calculate vital rate elasticieites #
#---------------------------------------#

## Using perturbation analysis
# TODO: 
# - Analogous to 6., but using formulas for elasticity (see 5.)

# Set perturbation factor
dy <- 1e-5

# Setting the perturbation for elasticity of target vital rates
E.pRep <- pRep * (1 + dy) # pertubate Breeding probability
E.mean.CS <- mean.CS * (1 + dy) # pertubate Mean clutch size
E.S_C <- S_C * (1 + dy) # pertubate Clutch survival
E.sF_NB_Ch <- sF_NB * c(1 + dy, 1, 1, 1) # pertubate Female chick Non-breeding survival
E.sF_NB_1yr <- sF_NB * c(1, 1 + dy, 1, 1) # pertubate Female 1yr Non-breeding survival
E.sF_NB_2yr <- sF_NB * c(1, 1, 1 + dy, 1) # pertubate Female 2yr Non-breeding survival
E.sF_NB_3yr <- sF_NB * c(1, 1, 1, 1 + dy) # pertubate Female 3yr Non-breeding survival
E.sF_BN_Ju <- sF_BN * c(1 + dy, 1, 1, 1) # pertubate Female Juvenile Breeding survival
E.sF_BN_1yr <- sF_BN * c(1, 1 + dy, 1, 1) # pertubate Female 1yr Breeding survival
E.sF_BN_2yr <- sF_BN * c(1, 1, 1 + dy, 1) # pertubate Female 2yr Breeding survival
E.sF_BN_3yr <- sF_BN * c(1, 1, 1, 1 + dy) # pertubate Female 3yr Breeding survival
E.sM_NB_Ch <- sM_NB * c(1 + dy, 1, 1, 1) # pertubate Male chick Non-breeding survival
E.sM_NB_1yr <- sM_NB * c(1, 1 + dy, 1, 1) # pertubate Male 1yr Non-breeding survival
E.sM_NB_2yr <- sM_NB * c(1, 1, 1 + dy, 1) # pertubate Male 2yr Non-breeding survival
E.sM_NB_3yr <- sM_NB * c(1, 1, 1, 1 + dy) # pertubate Male 3yr Non-breeding survival
E.sM_BN_Ju <- sM_BN * c(1 + dy, 1, 1, 1) # pertubate Male juvenile Breeding survival
E.sM_BN_1yr <- sM_BN * c(1, 1 + dy, 1, 1) # pertubate Male 1yr Breeding survival
E.sM_BN_2yr <- sM_BN * c(1, 1, 1 + dy, 1) # pertubate Male 2yr Breeding survival
E.sM_BN_3yr <- sM_BN * c(1, 1, 1, 1 + dy) # pertubate Male 3yr Breeding survival

# Origin Matrix
Ela_orig <- mat.ann$A

# Build matrix with pertubated target vital rate for elasticity

# Breeding probability
mat.Ela.pRep <- make.GPprojMatrix(pRep = E.pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                  sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                  seasonal = FALSE)

# Mean clutch size
mat.Ela.mean.CS <- make.GPprojMatrix(pRep = pRep, mean.CS = E.mean.CS, S_C = S_C, gamma = gamma, 
                                     sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                     seasonal = FALSE)

# Clutch survival
mat.Ela.S_C <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = E.S_C, gamma = gamma, 
                                 sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                 seasonal = FALSE)

# Female chick Non-breeding survival
mat.Ela.sF_NB_Ch <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = E.sF_NB_Ch, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                      seasonal = FALSE)

# Female 1yr Non-breeding survival
mat.Ela.sF_NB_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = E.sF_NB_1yr, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 2yr Non-breeding survival
mat.Ela.sF_NB_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = E.sF_NB_2yr, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 3yr Non-breeding survival
mat.Ela.sF_NB_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = E.sF_NB_3yr, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female Juvenile Breeding survival
mat.Ela.sF_BN_Ju <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = sF_NB, sF_BN = E.sF_BN_Ju, sM_NB = sM_NB, sM_BN = sM_BN,
                                      seasonal = FALSE)

# Female 1yr Breeding survival
mat.Ela.sF_BN_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = E.sF_BN_1yr, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 2yr Breeding survival
mat.Ela.sF_BN_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = E.sF_BN_2yr, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Female 3yr Breeding survival
mat.Ela.sF_BN_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = E.sF_BN_3yr, sM_NB = sM_NB, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male chick Non-breeding survival
mat.Ela.sM_NB_Ch <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = E.sM_NB_Ch, sM_BN = sM_BN,
                                      seasonal = FALSE)

# Male 1yr Non-breeding survival
mat.Ela.sM_NB_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = E.sM_NB_1yr, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male 2yr Non-breeding survival
mat.Ela.sM_NB_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = E.sM_NB_2yr, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male 3yr Non-breeding survival
mat.Ela.sM_NB_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = E.sM_NB_3yr, sM_BN = sM_BN,
                                       seasonal = FALSE)

# Male juvenile Breeding survival
mat.Ela.sM_BN_Ju <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                      sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = E.sM_BN_Ju,
                                      seasonal = FALSE)

# Male 1yr Breeding survival
mat.Ela.sM_BN_1yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = E.sM_BN_1yr,
                                       seasonal = FALSE)

# Male 2yr Breeding survival
mat.Ela.sM_BN_2yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = E.sM_BN_2yr,
                                       seasonal = FALSE)

# Male 3yr Breeding survival
mat.Ela.sM_BN_3yr <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                                       sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = E.sM_BN_3yr,
                                       seasonal = FALSE)

# Select only numerical matrix for calculation
Ela.pRep <- mat.Ela.pRep$A
Ela.mean.CS <- mat.Ela.mean.CS$A
Ela.S_C <- mat.Ela.S_C$A
Ela.sF_NB_Ch <- mat.Ela.sF_NB_Ch$A
Ela.sF_NB_1yr <- mat.Ela.sF_NB_1yr$A
Ela.sF_NB_2yr <- mat.Ela.sF_NB_2yr$A
Ela.sF_NB_3yr <- mat.Ela.sF_NB_3yr$A
Ela.sF_BN_Ju <- mat.Ela.sF_BN_Ju$A
Ela.sF_BN_1yr <- mat.Ela.sF_BN_1yr$A
Ela.sF_BN_2yr <- mat.Ela.sF_BN_2yr$A
Ela.sF_BN_3yr <- mat.Ela.sF_BN_3yr$A
Ela.sM_NB_Ch <- mat.Ela.sM_NB_Ch$A
Ela.sM_NB_1yr <- mat.Ela.sM_NB_1yr$A
Ela.sM_NB_2yr <- mat.Ela.sM_NB_2yr$A
Ela.sM_NB_3yr <- mat.Ela.sM_NB_3yr$A
Ela.sM_BN_Ju <- mat.Ela.sM_BN_Ju$A
Ela.sM_BN_1yr <- mat.Ela.sM_BN_1yr$A
Ela.sM_BN_2yr <- mat.Ela.sM_BN_2yr$A
Ela.sM_BN_3yr <- mat.Ela.sM_BN_3yr$A

# Calculate population growth rate for Origin matrix
lam_Ela_orig <- as.numeric(eigen(Ela_orig)$values[1])

# Calculate population growth rate for Origin matrix
lam_Ela.pRep <- as.numeric(eigen(Ela.pRep)$values[1])
lam_Ela.mean.CS <- as.numeric(eigen(Ela.mean.CS)$values[1])
lam_Ela.S_C <- as.numeric(eigen(Ela.S_C)$values[1])
lam_Ela.sF_NB_Ch <- as.numeric(eigen(Ela.sF_NB_Ch)$values[1])
lam_Ela.sF_NB_1yr <- as.numeric(eigen(Ela.sF_NB_1yr)$values[1])
lam_Ela.sF_NB_2yr <- as.numeric(eigen(Ela.sF_NB_2yr)$values[1])
lam_Ela.sF_NB_3yr <- as.numeric(eigen(Ela.sF_NB_3yr)$values[1])
lam_Ela.sF_BN_Ju <- as.numeric(eigen(Ela.sF_BN_Ju)$values[1])
lam_Ela.sF_BN_1yr <- as.numeric(eigen(Ela.sF_BN_1yr)$values[1])
lam_Ela.sF_BN_2yr <- as.numeric(eigen(Ela.sF_BN_2yr)$values[1])
lam_Ela.sF_BN_3yr <- as.numeric(eigen(Ela.sF_BN_3yr)$values[1])
lam_Ela.sM_NB_Ch <- as.numeric(eigen(Ela.sM_NB_Ch)$values[1])
lam_Ela.sM_NB_1yr <- as.numeric(eigen(Ela.sM_NB_1yr)$values[1])
lam_Ela.sM_NB_2yr <- as.numeric(eigen(Ela.sM_NB_2yr)$values[1])
lam_Ela.sM_NB_3yr <- as.numeric(eigen(Ela.sM_NB_3yr)$values[1])
lam_Ela.sM_BN_Ju <- as.numeric(eigen(Ela.sM_BN_Ju)$values[1])
lam_Ela.sM_BN_1yr <- as.numeric(eigen(Ela.sM_BN_1yr)$values[1])
lam_Ela.sM_BN_2yr <- as.numeric(eigen(Ela.sM_BN_2yr)$values[1])
lam_Ela.sM_BN_3yr <- as.numeric(eigen(Ela.sM_BN_3yr)$values[1])

# Calculate elasticity of population growth rate to target element
Ela_pRep <- (lam_Ela.pRep - lam_Ela_orig) / (lam_Ela_orig * dy) 
Ela_mean.CS <- (lam_Ela.mean.CS - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_S_C <- (lam_Ela.S_C - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_NB_Ch <- (lam_Ela.sF_NB_Ch - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_NB_1yr <- (lam_Ela.sF_NB_1yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_NB_2yr <- (lam_Ela.sF_NB_2yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_NB_3yr <- (lam_Ela.sF_NB_3yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_BN_Ju <- (lam_Ela.sF_BN_Ju - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_BN_1yr <- (lam_Ela.sF_BN_1yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_BN_2yr <- (lam_Ela.sF_BN_2yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sF_BN_3yr <- (lam_Ela.sF_BN_3yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_NB_Ch  <- (lam_Ela.sM_NB_Ch - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_NB_1yr <- (lam_Ela.sM_NB_1yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_NB_2yr <- (lam_Ela.sM_NB_2yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_NB_3yr <- (lam_Ela.sM_NB_3yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_BN_Ju <- (lam_Ela.sM_BN_Ju - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_BN_1yr <- (lam_Ela.sM_BN_1yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_BN_2yr <- (lam_Ela.sM_BN_2yr - lam_Ela_orig) / (lam_Ela_orig * dy)
Ela_sM_BN_3yr <- (lam_Ela.sM_BN_3yr - lam_Ela_orig) / (lam_Ela_orig * dy)

#Gather the results
Elasticities <- list(c(Ela_pRep, Ela_mean.CS, Ela_S_C, 
                       Ela_sF_NB_Ch, Ela_sF_NB_1yr, Ela_sF_NB_2yr, Ela_sF_NB_3yr,
                       Ela_sF_BN_Ju, Ela_sF_BN_1yr, Ela_sF_BN_2yr, Ela_sF_BN_3yr,
                       Ela_sM_NB_Ch, Ela_sM_NB_1yr, Ela_sM_NB_2yr, Ela_sM_NB_3yr,
                       Ela_sM_BN_Ju, Ela_sM_BN_1yr, Ela_sM_BN_2yr, Ela_sM_BN_3yr))




# - Visualize elasticities for all vital rates 






