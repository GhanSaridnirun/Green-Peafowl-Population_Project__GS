########################
# GREEN PEAFOWL IPM    #
# SENSITIVITY ANALYSIS #
########################

library(coda)

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


# Summarise sensitivity/elasticity arrays
# -> Median
# -> Mean
# -> Standard deviation

?apply()


