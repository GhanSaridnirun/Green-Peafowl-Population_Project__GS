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

A_pert[3, 1] <- A_pert[3, 1] + dy
A_pert[4, 2] <- A_pert[4, 2] + dy
A_pert[5, 3] <- A_pert[5, 3] + dy
A_pert[6, 4] <- A_pert[6, 4] + dy
A_pert[1, 5] <- A_pert[1, 5] + dy
A_pert[2, 5] <- A_pert[2, 5] + dy
A_pert[7, 5] <- A_pert[7, 5] + dy
A_pert[8, 6] <- A_pert[8, 6] + dy
A_pert[1, 7] <- A_pert[1, 7] + dy
A_pert[2, 7] <- A_pert[2, 7] + dy
A_pert[7, 7] <- A_pert[7, 7] + dy
A_pert[8, 8] <- A_pert[8, 8] + dy

A_pert

# Calculate population growth rate for both matrices
lam_orig <- as.numeric(eigen(A_orig)$values[1])
lam_pert <- as.numeric(eigen(A_pert)$values[1])

# Calculate sensitivity of population growth rate to target element
sens <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens[3, 1] <- (lam_pert - lam_orig) / dy 
sens[4, 2] <- (lam_pert - lam_orig) / dy 
sens[5, 3] <- (lam_pert - lam_orig) / dy 
sens[6, 4] <- (lam_pert - lam_orig) / dy 
sens[1, 5] <- (lam_pert - lam_orig) / dy 
sens[2, 5] <- (lam_pert - lam_orig) / dy 
sens[7, 5] <- (lam_pert - lam_orig) / dy 
sens[8, 6] <- (lam_pert - lam_orig) / dy 
sens[1, 7] <- (lam_pert - lam_orig) / dy 
sens[2, 7] <- (lam_pert - lam_orig) / dy 
sens[7, 7] <- (lam_pert - lam_orig) / dy 
sens[8, 8] <- (lam_pert - lam_orig) / dy 

sens


# TODO: 
# - Generalize above code to calculate sensitivities for all (non-0) matrix elements
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
A_pert[3, 1] <- A_pert[3, 1] * (1 + dy) 
A_pert[4, 2] <- A_pert[4, 2] * (1 + dy) 
A_pert[5, 3] <- A_pert[5, 3] * (1 + dy) 
A_pert[6, 4] <- A_pert[6, 4] * (1 + dy) 
A_pert[1, 5] <- A_pert[1, 5] * (1 + dy) 
A_pert[2, 5] <- A_pert[2, 5] * (1 + dy) 
A_pert[7, 5] <- A_pert[7, 5] * (1 + dy) 
A_pert[8, 6] <- A_pert[8, 6] * (1 + dy) 
A_pert[1, 7] <- A_pert[1, 7] * (1 + dy) 
A_pert[2, 7] <- A_pert[2, 7] * (1 + dy) 
A_pert[7, 7] <- A_pert[7, 7] * (1 + dy) 
A_pert[8, 8] <- A_pert[8, 8] * (1 + dy)

A_pert


# Calculate population growth rate for both matrices
lam_orig <- as.numeric(eigen(A_orig)$values[1])
lam_pert <- as.numeric(eigen(A_pert)$values[1])

# Calculate sensitivity of population growth rate to target element
elas <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas[3, 1] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[4, 2] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[5, 3] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[6, 4] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[1, 5] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[2, 5] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[7, 5] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[8, 6] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[1, 7] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[2, 7] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[7, 7] <- (lam_pert - lam_orig) / (lam_orig * dy) 
elas[8, 8] <- (lam_pert - lam_orig) / (lam_orig * dy) 

elas


# TODO: 
# - Generalize above code to calculate sensitivities for all (non-0) matrix elements
# - Visualize elasticity matrix (tip: start with fields::image.plot)

fields::image.plot(elas)



# 6. Calculate vital rate sensitivities #
#---------------------------------------#

# TODO: 
# - Adapt the above code to calculate sensitivities with respect to vital rates (instead of matrix elements)
#   NOTE: Keep in mind that unlike matrix elements, vital rates may appear several times in the matrix!
# - Generalize your code to calculate sensitivities for all vital rates

## Using perturbation analysis

# Origin Matrix
A_orig

# Perturbation Matrix
A_pert[3, 1] <- A_pert[3, 1] + dy 
A_pert[4, 2] <- A_pert[4, 2] + dy 
A_pert[5, 3] <- A_pert[5, 3] + dy 
A_pert[6, 4] <- A_pert[6, 4] + dy 
A_pert[1, 5] <- A_pert[1, 5] + dy 
A_pert[2, 5] <- A_pert[2, 5] + dy 
A_pert[7, 5] <- A_pert[7, 5] + dy 
A_pert[8, 6] <- A_pert[8, 6] + dy 
A_pert[1, 7] <- A_pert[1, 7] + dy 
A_pert[2, 7] <- A_pert[2, 7] + dy 
A_pert[7, 7] <- A_pert[7, 7] + dy 
A_pert[8, 8] <- A_pert[8, 8] + dy

A_pert 


# Calculate sensitivities with respect to vital rates

# Sensitivity with respect to Female chick survival  
Res.sF.Chick.orig <- A_orig/(sF_NB[1] * sF_BN[1])  
Res.sF.Chick.pert <- A_pert/(sF_NB[1] * sF_BN[1])

lam_Res.sF.Chick.orig <- as.numeric(eigen(Res.sF.Chick.orig)$values[1])
lam_Res.sF.Chick.pert <- as.numeric(eigen(Res.sF.Chick.pert)$values[1])

sens.Res.sF.Chick <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sF.Chick[3, 1] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[4, 2] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[5, 3] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[6, 4] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[1, 5] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[2, 5] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[7, 5] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[8, 6] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[1, 7] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[2, 7] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[7, 7] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy
sens.Res.sF.Chick[8, 8] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / dy

sens.Res.sF.Chick



# Sensitivity with respect to Male chick survival
Res.sM.Chick.orig <- A_orig/(sM_NB[1] * sM_BN[1]) 
Res.sM.Chick.pert <- A_pert/(sM_NB[1] * sM_BN[1])

lam_Res.sM.Chick.orig <- as.numeric(eigen(Res.sM.Chick.orig)$values[1])
lam_Res.sM.Chick.pert <- as.numeric(eigen(Res.sM.Chick.pert)$values[1])

sens.Res.sM.Chick <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sM.Chick[3, 1] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy 
sens.Res.sM.Chick[4, 2] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[5, 3] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy 
sens.Res.sM.Chick[6, 4] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy 
sens.Res.sM.Chick[1, 5] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[2, 5] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[7, 5] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[8, 6] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[1, 7] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[2, 7] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[7, 7] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy
sens.Res.sM.Chick[8, 8] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / dy

sens.Res.sM.Chick


# Sensitivity with respect to Female 1Y survival
Res.sF.1Y.orig <- A_orig/(sF_NB[2] * sF_BN[2])
Res.sF.1Y.pert <- A_pert/(sF_NB[2] * sF_BN[2])

lam_Res.sF.1Y.orig <- as.numeric(eigen(Res.sF.1Y.orig)$values[1])
lam_Res.sF.1Y.pert <- as.numeric(eigen(Res.sF.1Y.pert)$values[1])

sens.Res.sF.1Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sF.1Y[3, 1] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy 
sens.Res.sF.1Y[4, 2] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[5, 3] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[6, 4] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[1, 5] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[2, 5] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[7, 5] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[8, 6] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[1, 7] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[2, 7] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy
sens.Res.sF.1Y[7, 7] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy 
sens.Res.sF.1Y[8, 8] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / dy

sens.Res.sF.1Y


# Sensitivity with respect to Male 1Y survival
Res.sM.1Y.orig <- A_orig/(sM_NB[2] * sM_BN[2])
Res.sM.1Y.pert <- A_pert/(sM_NB[2] * sM_BN[2])

lam_Res.sM.1Y.orig <- as.numeric(eigen(Res.sM.1Y.orig)$values[1])
lam_Res.sM.1Y.pert <- as.numeric(eigen(Res.sM.1Y.pert)$values[1])

sens.Res.sM.1Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sM.1Y[3, 1] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[4, 2] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[5, 3] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[6, 4] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[1, 5] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[2, 5] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[7, 5] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[8, 6] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[1, 7] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[2, 7] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 
sens.Res.sM.1Y[7, 7] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy  
sens.Res.sM.1Y[8, 8] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / dy 

sens.Res.sM.1Y


# Sensitivity with respect to Female 2Y survival
Res.sF.2Y.orig <- A_orig/(sF_NB[3] * sF_BN[3])
Res.sF.2Y.pert <- A_pert/(sF_NB[3] * sF_BN[3])

lam_Res.sF.2Y.orig <- as.numeric(eigen(Res.sF.2Y.orig)$values[1])
lam_Res.sF.2Y.pert <- as.numeric(eigen(Res.sF.2Y.pert)$values[1])

sens.Res.sF.2Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sF.2Y[3, 1] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[4, 2] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[5, 3] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[6, 4] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[1, 5] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[2, 5] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[7, 5] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[8, 6] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[1, 7] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[2, 7] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 
sens.Res.sF.2Y[7, 7] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy  
sens.Res.sF.2Y[8, 8] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / dy 

sens.Res.sF.2Y


# Sensitivity with respect to Male 2Y survival
Res.sM.2Y.orig <- A_orig/(sM_NB[3] * sM_BN[3])
Res.sM.2Y.pert <- A_pert/(sM_NB[3] * sM_BN[3])

lam_Res.sM.2Y.orig <- as.numeric(eigen(Res.sM.2Y.orig)$values[1])
lam_Res.sM.2Y.pert <- as.numeric(eigen(Res.sM.2Y.pert)$values[1])

sens.Res.sM.2Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sM.2Y[3, 1] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[4, 2] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[5, 3] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[6, 4] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[1, 5] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[2, 5] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[7, 5] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[8, 6] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[1, 7] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 
sens.Res.sM.2Y[2, 7] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy  
sens.Res.sM.2Y[7, 7] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy  
sens.Res.sM.2Y[8, 8] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / dy 

sens.Res.sM.2Y


# Sensitivity with respect to Female 3+Y survival
Res.sF.3Y.orig <- A_orig/(sF_NB[4] * sF_BN[4])
Res.sF.3Y.pert <- A_pert/(sF_NB[4] * sF_BN[4])

lam_Res.sF.3Y.orig <- as.numeric(eigen(Res.sF.3Y.orig)$values[1])
lam_Res.sF.3Y.pert <- as.numeric(eigen(Res.sF.3Y.pert)$values[1])

sens.Res.sF.3Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sF.3Y[3, 1] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[4, 2] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[5, 3] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[6, 4] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[1, 5] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[2, 5] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[7, 5] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[8, 6] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[1, 7] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy 
sens.Res.sF.3Y[2, 7] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy  
sens.Res.sF.3Y[7, 7] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy  
sens.Res.sF.3Y[8, 8] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / dy  

sens.Res.sF.3Y


# Sensitivity with respect to Male 3+Y survival
Res.sM.3Y.orig <- A_orig/(sM_NB[4] * sM_BN[4])
Res.sM.3Y.pert <- A_pert/(sM_NB[4] * sM_BN[4])

lam_Res.sM.3Y.orig <- as.numeric(eigen(Res.sM.3Y.orig)$values[1])
lam_Res.sM.3Y.pert <- as.numeric(eigen(Res.sM.3Y.pert)$values[1])

sens.Res.sM.3Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.sM.3Y[3, 1] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[4, 2] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[5, 3] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[6, 4] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[1, 5] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy  
sens.Res.sM.3Y[2, 5] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[7, 5] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[8, 6] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[1, 7] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[2, 7] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy  
sens.Res.sM.3Y[7, 7] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy 
sens.Res.sM.3Y[8, 8] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / dy  

sens.Res.sM.3Y


# Sensitivity with respect to Breeding probability
Res.pRep.orig <- A_orig/pRep
Res.pRep.pert <- A_pert/pRep

lam_Res.pRep.orig <- as.numeric(eigen(Res.pRep.orig)$values[1])
lam_Res.pRep.pert <- as.numeric(eigen(Res.pRep.pert)$values[1])

sens.Res.pRep <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.pRep[3, 1] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy 
sens.Res.pRep[4, 2] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy  
sens.Res.pRep[5, 3] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy 
sens.Res.pRep[6, 4] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy 
sens.Res.pRep[1, 5] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy  
sens.Res.pRep[2, 5] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy  
sens.Res.pRep[7, 5] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy  
sens.Res.pRep[8, 6] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy 
sens.Res.pRep[1, 7] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy 
sens.Res.pRep[2, 7] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy  
sens.Res.pRep[7, 7] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy 
sens.Res.pRep[8, 8] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / dy  

sens.Res.pRep


# Sensitivity with respect to Mean clutch size
Res.mean.CS.orig <- A_orig/mean.CS
Res.mean.CS.pert <- A_pert/mean.CS

lam_Res.mean.CS.orig <- as.numeric(eigen(Res.mean.CS.orig)$values[1])
lam_Res.mean.CS.pert <- as.numeric(eigen(Res.mean.CS.pert)$values[1])

sens.Res.mean.CS <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.mean.CS[3, 1] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy 
sens.Res.mean.CS[4, 2] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy  
sens.Res.mean.CS[5, 3] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy 
sens.Res.mean.CS[6, 4] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy 
sens.Res.mean.CS[1, 5] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy  
sens.Res.mean.CS[2, 5] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy   
sens.Res.mean.CS[7, 5] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy  
sens.Res.mean.CS[8, 6] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy 
sens.Res.mean.CS[1, 7] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy  
sens.Res.mean.CS[2, 7] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy  
sens.Res.mean.CS[7, 7] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy 
sens.Res.mean.CS[8, 8] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / dy  

sens.Res.mean.CS


# Sensitivity with respect to Clutch survival
Res.S_C.orig <- A_orig/S_C
Res.S_C.pert <- A_pert/S_C

lam_Res.S_C.orig <- as.numeric(eigen(Res.S_C.orig)$values[1])
lam_Res.S_C.pert <- as.numeric(eigen(Res.S_C.pert)$values[1])

sens.Res.S_C <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
sens.Res.S_C[3, 1] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[4, 2] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[5, 3] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[6, 4] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy  
sens.Res.S_C[1, 5] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy   
sens.Res.S_C[2, 5] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy   
sens.Res.S_C[7, 5] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[8, 6] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[1, 7] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[2, 7] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[7, 7] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 
sens.Res.S_C[8, 8] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / dy 

sens.Res.S_C



# - Visualize sensitivities for all vital rates 

sense.all <- list(sens.Res.sF.Chick, sens.Res.sM.Chick, sens.Res.sF.1Y, sens.Res.sM.1Y,
                  sens.Res.sF.2Y, sens.Res.sM.2Y, sens.Res.sF.3Y, sens.Res.sM.1Y,
                  sens.Res.pRep, sens.Res.mean.CS, sens.Res.S_C)


# 7. Calculate vital rate elasticieites #
#---------------------------------------#

## Using perturbation analysis
# TODO: 
# - Analogous to 6., but using formulas for elasticity (see 5.)

# Origin Matrix
A_orig.elas <- A_orig

# Perturbation Matrix
A_pert.elas <- A_pert

A_pert.elas[3, 1] <- A_pert.elas[3, 1] * (1 + dy) # ChF to 1yF
A_pert.elas[4, 2] <- A_pert.elas[4, 2] * (1 + dy) # ChM to 1yM
A_pert.elas[5, 3] <- A_pert.elas[5, 3] * (1 + dy)# 1yF to 2yF
A_pert.elas[6, 4] <- A_pert.elas[6, 4] * (1 + dy) # 1yM to 2yM
A_pert.elas[1, 5] <- A_pert.elas[1, 5] * (1 + dy) # 2yF produce ChF
A_pert.elas[2, 5] <- A_pert.elas[2, 5] * (1 + dy) # 2yF produce ChM
A_pert.elas[7, 5] <- A_pert.elas[7, 5] * (1 + dy) # 2yF to 3yF
A_pert.elas[8, 6] <- A_pert.elas[8, 6] * (1 + dy)# 2yM to 3yM
A_pert.elas[1, 7] <- A_pert.elas[1, 7] * (1 + dy) # 3yF produce ChF
A_pert.elas[2, 7] <- A_pert.elas[2, 7] * (1 + dy) # 3yF produce ChM
A_pert.elas[7, 7] <- A_pert.elas[7, 7] * (1 + dy)# 3yF to 3yF
A_pert.elas[8, 8] <- A_pert.elas[8, 8] * (1 + dy)# 3yM to 3yM

A_pert.elas


# Calculate elasticities with respect to vital rates

# Elasticity with respect to Female chick survival  
Res.elas.sF.Chick.orig <- A_orig.elas/(sF_NB[1] * sF_BN[1])  
Res.elas.sF.Chick.pert <- A_pert.elas/(sF_NB[1] * sF_BN[1])

lam_Res.sF.Chick.orig <- as.numeric(eigen(Res.sF.Chick.orig)$values[1])
lam_Res.sF.Chick.pert <- as.numeric(eigen(Res.sF.Chick.pert)$values[1])

elas.Res.sF.Chick <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sF.Chick[3, 1] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # ChF to 1yF
elas.Res.sF.Chick[4, 2] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # ChM to 1yM
elas.Res.sF.Chick[5, 3] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 1yF to 2yF
elas.Res.sF.Chick[6, 4] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 1yM to 2yM
elas.Res.sF.Chick[1, 5] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 2yF produce ChF
elas.Res.sF.Chick[2, 5] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 2yF produce ChM
elas.Res.sF.Chick[7, 5] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 2yF to 3yF
elas.Res.sF.Chick[8, 6] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 2yM to 3yM
elas.Res.sF.Chick[1, 7] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 3yF produce ChF
elas.Res.sF.Chick[2, 7] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 3yF produce ChM
elas.Res.sF.Chick[7, 7] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 3yF to 3yF
elas.Res.sF.Chick[8, 8] <- (lam_Res.sF.Chick.pert - lam_Res.sF.Chick.orig) / (lam_Res.sF.Chick.orig * dy) # 3yM to 3yM

elas.Res.sF.Chick


# Sensitivity with respect to Male chick survival
Res.elas.sM.Chick.orig <- A_orig.elas/(sM_NB[1] * sM_BN[1])  
Res.elas.sM.Chick.pert <- A_pert.elas/(sM_NB[1] * sM_BN[1])

lam_Res.sM.Chick.orig <- as.numeric(eigen(Res.sM.Chick.orig)$values[1])
lam_Res.sM.Chick.pert <- as.numeric(eigen(Res.sM.Chick.pert)$values[1])

elas.Res.sM.Chick <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sM.Chick[3, 1] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[4, 2] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[5, 3] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[6, 4] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[1, 5] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[2, 5] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[7, 5] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[8, 6] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[1, 7] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[2, 7] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[7, 7] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)
elas.Res.sM.Chick[8, 8] <- (lam_Res.sM.Chick.pert - lam_Res.sM.Chick.orig) / (lam_Res.sM.Chick.orig * dy)

elas.Res.sM.Chick


# Sensitivity with respect to Female 1Y survival
Res.elas.sF.1Y.orig <- A_orig.elas/(sF_NB[2] * sF_BN[2])  
Res.elas.sF.1Y.pert <- A_pert.elas/(sF_NB[2] * sF_BN[2])

lam_Res.sF.1Y.orig <- as.numeric(eigen(Res.elas.sF.1Y.orig )$values[1])
lam_Res.sF.1Y.pert <- as.numeric(eigen(Res.elas.sF.1Y.pert)$values[1])

elas.Res.sF.1Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sF.1Y[3, 1] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[4, 2] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[5, 3] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[6, 4] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[1, 5] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[2, 5] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[7, 5] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[8, 6] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[1, 7] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[2, 7] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[7, 7] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)
elas.Res.sF.1Y[8, 8] <- (lam_Res.sF.1Y.pert - lam_Res.sF.1Y.orig) / (lam_Res.sF.1Y.orig * dy)

elas.Res.sF.1Y


# Sensitivity with respect to Male 1Y survival
Res.elas.sM.1Y.orig <- A_orig.elas/(sM_NB[2] * sM_BN[2])  
Res.elas.sM.1Y.pert <- A_pert.elas/(sM_NB[2] * sM_BN[2])

lam_Res.sM.1Y.orig <- as.numeric(eigen(Res.elas.sM.1Y.orig )$values[1])
lam_Res.sM.1Y.pert <- as.numeric(eigen(Res.elas.sM.1Y.pert)$values[1])

elas.Res.sM.1Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sM.1Y[3, 1] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[4, 2] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[5, 3] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[6, 4] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[1, 5] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[2, 5] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[7, 5] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[8, 6] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[1, 7] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[2, 7] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[7, 7] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)
elas.Res.sM.1Y[8, 8] <- (lam_Res.sM.1Y.pert - lam_Res.sM.1Y.orig) / (lam_Res.sM.1Y.orig * dy)

elas.Res.sM.1Y


# Sensitivity with respect to Female 2Y survival
Res.elas.sF.2Y.orig <- A_orig.elas/(sF_NB[3] * sF_BN[3])  
Res.elas.sF.2Y.pert <- A_pert.elas/(sF_NB[3] * sF_BN[3])

lam_Res.sF.2Y.orig <- as.numeric(eigen(Res.elas.sF.2Y.orig )$values[1])
lam_Res.sF.2Y.pert <- as.numeric(eigen(Res.elas.sF.2Y.pert)$values[1])

elas.Res.sF.2Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sF.2Y[3, 1] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[4, 2] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[5, 3] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[6, 4] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[1, 5] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[2, 5] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[7, 5] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[8, 6] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[1, 7] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[2, 7] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[7, 7] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)
elas.Res.sF.2Y[8, 8] <- (lam_Res.sF.2Y.pert - lam_Res.sF.2Y.orig) / (lam_Res.sF.2Y.orig * dy)

elas.Res.sF.2Y


# Sensitivity with respect to Male 2Y survival
Res.elas.sM.2Y.orig <- A_orig.elas/(sM_NB[3] * sM_BN[3])  
Res.elas.sM.2Y.pert <- A_pert.elas/(sM_NB[3] * sM_BN[3])

lam_Res.sM.2Y.orig <- as.numeric(eigen(Res.elas.sM.2Y.orig )$values[1])
lam_Res.sM.2Y.pert <- as.numeric(eigen(Res.elas.sM.2Y.pert)$values[1])

elas.Res.sM.2Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sM.2Y[3, 1] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[4, 2] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[5, 3] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[6, 4] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[1, 5] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[2, 5] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[7, 5] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[8, 6] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[1, 7] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[2, 7] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[7, 7] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)
elas.Res.sM.2Y[8, 8] <- (lam_Res.sM.2Y.pert - lam_Res.sM.2Y.orig) / (lam_Res.sM.2Y.orig * dy)

elas.Res.sM.2Y


# Sensitivity with respect to Female 3+Y survival
Res.elas.sF.3Y.orig <- A_orig.elas/(sF_NB[4] * sF_BN[4])  
Res.elas.sF.3Y.pert <- A_pert.elas/(sF_NB[4] * sF_BN[4])

lam_Res.sF.3Y.orig <- as.numeric(eigen(Res.elas.sF.3Y.orig )$values[1])
lam_Res.sF.3Y.pert <- as.numeric(eigen(Res.elas.sF.3Y.pert)$values[1])

elas.Res.sF.3Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sF.3Y[3, 1] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[4, 2] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[5, 3] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[6, 4] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[1, 5] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[2, 5] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[7, 5] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[8, 6] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[1, 7] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[2, 7] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[7, 7] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)
elas.Res.sF.3Y[8, 8] <- (lam_Res.sF.3Y.pert - lam_Res.sF.3Y.orig) / (lam_Res.sF.3Y.orig * dy)

elas.Res.sF.3Y


# Sensitivity with respect to Male 3+Y survival
Res.elas.sM.3Y.orig <- A_orig.elas/(sM_NB[4] * sM_BN[4])  
Res.elas.sM.3Y.pert <- A_pert.elas/(sM_NB[4] * sM_BN[4])

lam_Res.sM.3Y.orig <- as.numeric(eigen(Res.elas.sM.3Y.orig )$values[1])
lam_Res.sM.3Y.pert <- as.numeric(eigen(Res.elas.sM.3Y.pert)$values[1])

elas.Res.sM.3Y <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.sM.3Y[3, 1] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[4, 2] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[5, 3] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[6, 4] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[1, 5] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[2, 5] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[7, 5] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[8, 6] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[1, 7] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[2, 7] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[7, 7] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)
elas.Res.sM.3Y[8, 8] <- (lam_Res.sM.3Y.pert - lam_Res.sM.3Y.orig) / (lam_Res.sM.3Y.orig * dy)

elas.Res.sM.3Y


# Sensitivity with respect to Breeding probability
Res.elas.pRep.orig <- A_orig.elas/pRep  
Res.elas.pRep.pert <- A_pert.elas/pRep

lam_Res.pRep.orig <- as.numeric(eigen(Res.elas.pRep.orig )$values[1])
lam_Res.pRep.pert <- as.numeric(eigen(Res.elas.pRep.pert)$values[1])

elas.Res.pRep <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.pRep[3, 1] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[4, 2] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[5, 3] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[6, 4] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[1, 5] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[2, 5] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[7, 5] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[8, 6] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[1, 7] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[2, 7] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[7, 7] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)
elas.Res.pRep[8, 8] <- (lam_Res.pRep.pert - lam_Res.pRep.orig) / (lam_Res.pRep.orig * dy)

elas.Res.pRep


# Sensitivity with respect to Mean clutch size
Res.elas.mean.CS.orig <- A_orig.elas/mean.CS
Res.elas.mean.CS.pert <- A_pert.elas/mean.CS

lam_Res.mean.CS.orig <- as.numeric(eigen(Res.elas.mean.CS.orig )$values[1])
lam_Res.mean.CS.pert <- as.numeric(eigen(Res.elas.mean.CS.pert)$values[1])

elas.Res.mean.CS <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.mean.CS[3, 1] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[4, 2] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[5, 3] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[6, 4] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[1, 5] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[2, 5] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[7, 5] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[8, 6] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[1, 7] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[2, 7] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[7, 7] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)
elas.Res.mean.CS[8, 8] <- (lam_Res.mean.CS.pert - lam_Res.mean.CS.orig) / (lam_Res.mean.CS.orig * dy)

elas.Res.mean.CS


# Sensitivity with respect to Clutch survival
Res.elas.S_C.orig <- A_orig.elas/S_C
Res.elas.S_C.pert <- A_pert.elas/S_C

lam_Res.S_C.orig <- as.numeric(eigen(Res.elas.S_C.orig )$values[1])
lam_Res.S_C.pert <- as.numeric(eigen(Res.elas.S_C.pert)$values[1])

elas.Res.S_C <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig))
elas.Res.S_C[3, 1] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[4, 2] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[5, 3] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[6, 4] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[1, 5] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[2, 5] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[7, 5] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[8, 6] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[1, 7] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[2, 7] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[7, 7] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)
elas.Res.S_C[8, 8] <- (lam_Res.S_C.pert - lam_Res.S_C.orig) / (lam_Res.S_C.orig * dy)

elas.Res.S_C

# - Visualize elasticities for all vital rates 

elas.all <- list(elas.Res.sF.Chick, elas.Res.sM.Chick, elas.Res.sF.1Y, elas.Res.sM.1Y,
                 elas.Res.sF.2Y, elas.Res.sM.2Y, elas.Res.sF.3Y, elas.Res.sM.1Y,
                 elas.Res.pRep, elas.Res.mean.CS, elas.Res.S_C)





