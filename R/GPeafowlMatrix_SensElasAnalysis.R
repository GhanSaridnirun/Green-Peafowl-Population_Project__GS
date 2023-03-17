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

# Set up original and sensitivity matrices
A_orig <- mat.ann$A
sens <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig), dimnames = dimnames(A_orig))

calculate.sensitivity <- function(A_orig, dy){

                for(j in 1:nrow(A_orig)){
                  for(k in 1:ncol(A_orig)){
    
                      if(A_orig[j, k] != 0){ # Only continue if matrix element [j, k] != 0
      
                        # Define perturbation matrix  
                        A_pert <- A_orig
      
                        # Perturb target element in matrix
                        A_pert[j, k] <- A_pert[j, k] + dy
      
                        # Calculate population growth rate for both matrices
                        lam_orig <- as.numeric(eigen(A_orig)$values[1])
                        lam_pert <- as.numeric(eigen(A_pert)$values[1])
      
                        # Calculate sensitivity of population growth rate to target element
                        sens[j, k] <- (lam_pert - lam_orig) / dy 
                      }
                    }
                  }

    return(sens)
} 

sens <- calculate.sensitivity(A_orig, dy)

library(fields)

## Visualize sensitivity matrix
fields::image.plot(t(apply(sens, 2, rev)), axes = FALSE, col = plasma(20))

image(t(apply(sens, 2, rev)), col = plasma(20), add = TRUE)
axis(3, at = seq(0, 1, length = ncol(sens)), labels = colnames(sens), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(sens)), labels = rownames(sens), lwd = 0, las = 2)


# 5. Calculate matrix element elasticities #
#------------------------------------------#

## Using perturbation analysis

# Set perturbation factor
dy <- 1e-5

# Set up original and perturbed matrix
A_orig <- mat.ann$A
elas <- matrix(0, nrow = nrow(A_orig), ncol = ncol(A_orig), dimnames = dimnames(A_orig))

calculate.elasticity <- function(A_orig, dy){

               for(m in 1:nrow(A_orig)){
                for(n in 1:ncol(A_orig)){
    
                    if(A_orig[m, n] != 0){ # Only continue if matrix element [m, n] != 0
      
                      # Define perturbation matrix  
                      A_pert <- A_orig
      
                      # Perturb target element in matrix
                      A_pert[m, n] <- A_pert[m, n] * (1 + dy) 
      
                      # Calculate population growth rate for both matrices
                      lam_orig <- as.numeric(eigen(A_orig)$values[1])
                      lam_pert <- as.numeric(eigen(A_pert)$values[1])
      
                      # Calculate sensitivity of population growth rate to target element
                      elas[m, n] <- (lam_pert - lam_orig) / (lam_orig * dy)
                    }
                  }
                }
  return(elas)
}  


elas <- calculate.elasticity(A_orig, dy)


# - Visualize elasticity matrix (tip: start with fields::image.plot)

## Visualize elasticity matrix
fields::image.plot(t(apply(elas, 2, rev)), axes = FALSE, col = plasma(20))

image(t(apply(elas, 2, rev)), col = plasma(20), add = TRUE)
axis(3, at = seq(0, 1, length = ncol(elas)), labels = colnames(elas), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(elas)), labels = rownames(elas), lwd = 0, las = 2)



# 6. Calculate vital rate sensitivities #
#---------------------------------------#

## Using perturbation analysis

# Set perturbation factor
dy <- 1e-5

# List all vital rate names
VR.names <- c("pRep", "mean.CS", "gamma", "S_C",
              "sF_NB_Ch", paste0("sF_NB_", 1:3, "yr"),
              "sF_BN_Ju", paste0("sF_BN_", 1:3, "yr"),
              "sM_NB_Ch", paste0("sM_NB_", 1:3, "yr"),
              "sM_BN_Ju", paste0("sM_BN_", 1:3, "yr"))

# List all original vital rate values
VR.orig <- c(pRep, mean.CS, gamma, S_C,
             sF_NB, sF_BN, sM_NB, sM_BN)
names(VR.orig) <- VR.names

# Make perturbation matrix
# NOTE: Each row in this matrix corresponds to one "perturbation scenario", 
# i.e. a situation in which the vital rate with the corresponding row name is 
# perturbed. By having each scenario in a column, we can loop over all the 
# possibilities instead of writing them out manually. 
pert.mat <- diag(length(VR.names))*dy
rownames(pert.mat) <- VR.names

# Make matrix of perturbed vital rates (per scenario, sensitivity = additive)
VR.pert <- VR.orig + pert.mat

# Build original matrix
A_orig <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                            sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                            seasonal = FALSE)$A

# Set up empty dataframe for storing results
VR.sens.data <- data.frame()

calculate.VR.Sens <- function(dy, VR.names, VR.orig){

                      for(i in 1:length(VR.names)){ # Loop over each perturbation scenario
  
                        # Build perturbed matrix
                        A_pert <- make.GPprojMatrix(pRep = VR.pert["pRep", i], 
                                                    mean.CS = VR.pert["mean.CS", i], 
                                                    S_C = VR.pert["S_C", i], 
                                                    gamma = VR.pert["gamma", i], 
                                                    sF_NB = VR.pert[c("sF_NB_Ch", paste0("sF_NB_", 1:3, "yr")), i], 
                                                    sF_BN = VR.pert[c("sF_BN_Ju", paste0("sF_BN_", 1:3, "yr")), i], 
                                                    sM_NB = VR.pert[c("sM_NB_Ch", paste0("sM_NB_", 1:3, "yr")), i], 
                                                    sM_BN = VR.pert[c("sM_BN_Ju", paste0("sM_BN_", 1:3, "yr")), i],
                                                    seasonal = FALSE)$A
  
                        # Calculate population growth rate for both matrices
                        lam_orig <- as.numeric(eigen(A_orig)$values[1])
                        lam_pert <- as.numeric(eigen(A_pert)$values[1])
  
                        # Calculate sensitivity of population growth rate to target element
                        VR.sens <- (lam_pert - lam_orig) / dy
  
                        # Assemble results in a dataframe
                        data.temp <- data.frame(VitalRate = VR.names[i],
                                                Sensitivity = VR.sens)
  
                        VR.sens.data <- rbind(VR.sens.data, data.temp)
                      }
  return(VR.sens.data)
}

VR.sens.data <- calculate.VR.Sens(dy, VR.names, VR.orig)

# Visualize sensitivities for all vital rates 
library(ggplot2)
ggplot(VR.sens.data) + 
  geom_bar(aes(x = VitalRate, y = Sensitivity), stat = "identity") + 
  coord_flip() + 
  theme_classic()

# TODO 1: Turn the above calculation of sensitivity (lines 266-312) into a 
# function that takes dy, VR.names, and VR.orig as input and returns VR.sens.data
# TODO 2: Write equivalent code for elasticities below


# 7. Calculate vital rate elasticiites #
#---------------------------------------#

## Using perturbation analysis
# - Analogous to 6., but using formulas for elasticity (see 5.)

# List all vital rate names
VR.names <- c("pRep", "mean.CS", "gamma", "S_C",
              "sF_NB_Ch", paste0("sF_NB_", 1:3, "yr"),
              "sF_BN_Ju", paste0("sF_BN_", 1:3, "yr"),
              "sM_NB_Ch", paste0("sM_NB_", 1:3, "yr"),
              "sM_BN_Ju", paste0("sM_BN_", 1:3, "yr"))

# List all original vital rate values
VR.orig <- c(pRep, mean.CS, gamma, S_C,
             sF_NB, sF_BN, sM_NB, sM_BN)
names(VR.orig) <- VR.names

# Make perturbation matrix
# NOTE: Each row in this matrix corresponds to one "perturbation scenario", 
# i.e. a situation in which the vital rate with the corresponding row name is 
# perturbed. By having each scenario in a column, we can loop over all the 
# possibilities instead of writing them out manually. 
pert.mat <- diag(length(VR.names))*dy
rownames(pert.mat) <- VR.names

# Make matrix of perturbed vital rates (per scenario, sensitivity = additive)
VR.pert <- VR.orig + pert.mat

# Build original matrix
A_orig <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                            sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                            seasonal = FALSE)$A

# Set up empty dataframe for storing results
VR.elas.data <- data.frame()

calculate.VR.Elas <- function(dy, VR.names, VR.orig){
  
                       for(i in 1:length(VR.names)){ # Loop over each perturbation scenario
    
                         # Build perturbed matrix
                         A_pert <- make.GPprojMatrix(pRep = VR.pert["pRep", i], 
                                                     mean.CS = VR.pert["mean.CS", i], 
                                                     S_C = VR.pert["S_C", i], 
                                                     gamma = VR.pert["gamma", i], 
                                                     sF_NB = VR.pert[c("sF_NB_Ch", paste0("sF_NB_", 1:3, "yr")), i], 
                                                     sF_BN = VR.pert[c("sF_BN_Ju", paste0("sF_BN_", 1:3, "yr")), i], 
                                                     sM_NB = VR.pert[c("sM_NB_Ch", paste0("sM_NB_", 1:3, "yr")), i], 
                                                     sM_BN = VR.pert[c("sM_BN_Ju", paste0("sM_BN_", 1:3, "yr")), i],
                                                     seasonal = FALSE)$A
    
                        # Calculate population growth rate for both matrices
                        lam_orig <- as.numeric(eigen(A_orig)$values[1])
                        lam_pert <- as.numeric(eigen(A_pert)$values[1])
    
                        # Calculate sensitivity of population growth rate to target element
                        VR.elas <- (lam_pert - lam_orig) / (lam_orig * dy)
    
                        # Assemble results in a dataframe
                        data.temp <- data.frame(VitalRate = VR.names[i],
                                                Elasticity = VR.elas)
    
                        VR.elas.data <- rbind(VR.elas.data, data.temp)
                      }
  return(VR.elas.data)
}

VR.elas.data <- calculate.VR.Elas(dy, VR.names, VR.orig)

# Visualize sensitivities for all vital rates 
library(ggplot2)
ggplot(VR.elas.data) + 
  geom_bar(aes(x = VitalRate, y = Elasticity), stat = "identity") + 
  coord_flip() + 
  theme_classic()










