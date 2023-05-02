########################
# GREEN PEAFOWL IPM    #
# SENSITIVITY ANALYSIS #
########################

library(coda)
library(viridis)

# Load functions for sensitivity and elasticity calculations
source("R/make.GPprojMatrix.R")
source("R/calculate.sensitivity.R")
source("R/calculate.elasticity.R")
source("R/calculate.VR.Sens.R")
source("R/calculate.VR.Elas.R")

# Read data from rds
IPM.rhoDeriv <- readRDS("GPIPM_TwoSex_Matrix_Clutch_BreedProb_rhoDeriv.rds")

# Set as matrix
mat.ipm <- as.matrix(IPM.rhoDeriv)

# Set number of samples
nSamples <- nrow(mat.ipm)

# Set perturbation factor
dy <- 1e-5

# Prepare objects for storing results
lambda.ann <- rep(NA, nSamples)
sens <- elas <- array(NA, dim = c(8, 8, nSamples))
VR.sens <- VR.elas <- data.frame()


for(i in 1:nSamples){
  
  
  # 1. Set parameter values #
  #-------------------------#
  
  ## Productivity 
  
  pRep <- mat.ipm[i, "pRep"]   
  mean.CS <- mat.ipm[i, "mean.CS"]
  gamma <- 0.5 
  S_C <- mat.ipm[i, "mean.S_C"]
  
  
  ## Survival
  
  sF_NB <- c(mat.ipm[i, "sF_NB[1]"], mat.ipm[i, "sF_NB[2]"], 
             mat.ipm[i, "sF_NB[3]"], mat.ipm[i, "sF_NB[4]"])
  
  sF_BN <- c(mat.ipm[i, "sF_BN[1]"], mat.ipm[i, "sF_BN[2]"], 
             mat.ipm[i, "sF_BN[3]"], mat.ipm[i, "sF_BN[4]"])
  
  sM_NB <- c(mat.ipm[i, "sM_NB[1]"], mat.ipm[i, "sM_NB[2]"], 
             mat.ipm[i, "sM_NB[3]"], mat.ipm[i, "sM_NB[4]"])
  
  sM_BN <- c(mat.ipm[i, "sM_BN[1]"], mat.ipm[i, "sM_BN[2]"], 
             mat.ipm[i, "sM_BN[3]"], mat.ipm[i, "sM_BN[4]"])
  
  
  
  # 2. Build projection matrix #
  #----------------------------#
  
  
  ## Build annual matrix
  mat.ann <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                               sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                               seasonal = FALSE)
  
  
  # 3. Calculate asymptotic population growth rate #
  #------------------------------------------------#
  
  
  ## Annual matrix
  # Using eigen decomposition
  lambda.ann[i] <- as.numeric(eigen(mat.ann$A)$values[1])
  
  
  
  # 4. Calculate matrix element sensitivities #
  #-------------------------------------------#
  
  
  # Set up original and sensitivity matrices
  A_orig <- mat.ann$A

  # Calculate Sensitivities
  sens[, , i] <- calculate.sensitivity(A_orig, dy)
  
  
  
  # 5. Calculate matrix element elasticities #
  #------------------------------------------#
  
  ## Using perturbation analysis
  
  # Set up original and perturbed matrix
  A_orig <- mat.ann$A

  # Calculate Elasticities
  elas[, , i] <- calculate.elasticity(A_orig, dy)
  
  
  
  # 6. Calculate vital rate sensitivities #
  #---------------------------------------#
  
  ## Using perturbation analysis
  
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
  
  # Calculate.VR.Sens 
  VR.sens.data <- calculate.VR.Sens(dy, VR.names, VR.orig)
  
  # Add sample marker
  VR.sens.data$Sample <- i
  
  # Bind to previous results
  VR.sens <- rbind(VR.sens, VR.sens.data)
  
  
  # 7. Calculate vital rate elasticieites #
  #---------------------------------------#
  
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
  rownames(pert.mat) <- VR.names
  
  # Make matrix of perturbed vital rates (per scenario, sensitivity = additive)
  VR.pert <- VR.orig + pert.mat
  
  # Build original matrix
  A_orig <- make.GPprojMatrix(pRep = pRep, mean.CS = mean.CS, S_C = S_C, gamma = gamma, 
                              sF_NB = sF_NB, sF_BN = sF_BN, sM_NB = sM_NB, sM_BN = sM_BN,
                              seasonal = FALSE)$A
  
  # Set up empty dataframe for storing results
  VR.elas.data <- data.frame()
  
  # Calculate.VR.Elas
  VR.elas.data <- calculate.VR.Elas(dy, VR.names, VR.orig)
  
  # Add sample marker
  VR.elas.data$Sample <- i
  
  # Bind to previous results
  VR.elas <- rbind(VR.elas, VR.elas.data)
  
}


## Summarise sensitivity/elasticity arrays

# Median
sens.median <- apply(sens, MARGIN = c(1, 2), FUN = median)
elas.median <- apply(elas, MARGIN = c(1, 2), FUN = median)

# Mean
sens.mean <- apply(sens, MARGIN = c(1, 2), FUN = mean)
elas.mean <- apply(elas, MARGIN = c(1, 2), FUN = mean)

# Standard deviation
sens.sd <- apply(sens, MARGIN = c(1, 2), FUN = sd)
elas.sd <- apply(elas, MARGIN = c(1, 2), FUN = sd)
 
## Heat plots of summaries

# Sensitivities - median
fields::image.plot(t(apply(sens.median, 2, rev)), axes = FALSE, col = plasma(20))
axis(3, at = seq(0, 1, length = ncol(A_orig)), labels = colnames(A_orig), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(A_orig)), labels = rownames(A_orig), lwd = 0, las = 2)

# Sensitivities - mean
fields::image.plot(t(apply(sens.mean, 2, rev)), axes = FALSE, col = plasma(20))
axis(3, at = seq(0, 1, length = ncol(A_orig)), labels = colnames(A_orig), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(A_orig)), labels = rownames(A_orig), lwd = 0, las = 2)

# Sensitivities - standard deviation (= uncertainty)
fields::image.plot(t(apply(sens.sd, 2, rev)), axes = FALSE, col = rev(grey.colors(20)))
axis(3, at = seq(0, 1, length = ncol(A_orig)), labels = colnames(A_orig), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(A_orig)), labels = rownames(A_orig), lwd = 0, las = 2)


# Elasticities - median
fields::image.plot(t(apply(elas.median, 2, rev)), axes = FALSE, col = plasma(20))
axis(3, at = seq(0, 1, length = ncol(A_orig)), labels = colnames(A_orig), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(A_orig)), labels = rownames(A_orig), lwd = 0, las = 2)

# Elasticities - mean
fields::image.plot(t(apply(elas.mean, 2, rev)), axes = FALSE, col = plasma(20))
axis(3, at = seq(0, 1, length = ncol(A_orig)), labels = colnames(A_orig), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(A_orig)), labels = rownames(A_orig), lwd = 0, las = 2)

# Elasticities - standard deviation (= uncertainty)
fields::image.plot(t(apply(elas.sd, 2, rev)), axes = FALSE, col = rev(grey.colors(20)))
axis(3, at = seq(0, 1, length = ncol(A_orig)), labels = colnames(A_orig), lwd = 0)
axis(2, at = seq(1, 0, length = nrow(A_orig)), labels = rownames(A_orig), lwd = 0, las = 2)




