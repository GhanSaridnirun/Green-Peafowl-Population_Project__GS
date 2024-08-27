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



# The extension part for process with sensitivity and elasticity 
# of vital rates

## Data in sensitivity and elasticity of Vital rates

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Select vital rates for calculate median, mean and sd of sensitivity

pRep.VR.S <- subset(VR.sens, VitalRate %in% "pRep")
mean.CS.VR.S <- subset(VR.sens, VitalRate %in% "mean.CS")
S_C.VR.S <- subset(VR.sens, VitalRate %in% "S_C")
sF_NB_Ch.VR.S <- subset(VR.sens, VitalRate %in% "sF_NB_Ch")
sF_NB_1yr.VR.S <- subset(VR.sens, VitalRate %in% "sF_NB_1yr")
sF_NB_2yr.VR.S <- subset(VR.sens, VitalRate %in% "sF_NB_2yr")
sF_NB_3yr.VR.S <- subset(VR.sens, VitalRate %in% "sF_NB_3yr")
sF_BN_Ju.VR.S <- subset(VR.sens, VitalRate %in% "sF_BN_Ju")
sF_BN_1yr.VR.S <- subset(VR.sens, VitalRate %in% "sF_BN_1yr")
sF_BN_2yr.VR.S <- subset(VR.sens, VitalRate %in% "sF_BN_2yr")
sF_BN_3yr.VR.S <- subset(VR.sens, VitalRate %in% "sF_BN_3yr")


# Select vital rates for calculate median, mean and sd of elasticity

pRep.VR.E <- subset(VR.elas, VitalRate %in% "pRep")
mean.CS.VR.E <- subset(VR.elas, VitalRate %in% "mean.CS")
S_C.VR.E <- subset(VR.elas, VitalRate %in% "S_C")
sF_NB_Ch.VR.E <- subset(VR.elas, VitalRate %in% "sF_NB_Ch")
sF_NB_1yr.VR.E <- subset(VR.elas, VitalRate %in% "sF_NB_1yr")
sF_NB_2yr.VR.E <- subset(VR.elas, VitalRate %in% "sF_NB_2yr")
sF_NB_3yr.VR.E <- subset(VR.elas, VitalRate %in% "sF_NB_3yr")
sF_BN_Ju.VR.E <- subset(VR.elas, VitalRate %in% "sF_BN_Ju")
sF_BN_1yr.VR.E <- subset(VR.elas, VitalRate %in% "sF_BN_1yr")
sF_BN_2yr.VR.E <- subset(VR.elas, VitalRate %in% "sF_BN_2yr")
sF_BN_3yr.VR.E <- subset(VR.elas, VitalRate %in% "sF_BN_3yr")



#Prepare data frame to store results.

sum.vr.sensitivity <- sum.vr.elasticity <- data.frame(
  Vitalrate = c("Breeding probability", "Clutch Size", "Clutch Survival",
                "Female Chick survival", 
                "Female Yearling Non-breeding Survival", 
                "Female 2 Years Non-breeding Survival", 
                "Female 3 Years Non-breeding Survival", 
                "Female Juvenile survival",
                "Female Yearling Breeding Survival", 
                "Female 2 Years Breeding Survival",
                "Female 3 Years Breeding Survival"),
  Median = NA,
  Mean = NA,
  SD = NA
)                      



# Calculate median, mean and sd of sensitivity of Vital rates

# Breeding probability
sum.vr.sensitivity[1,2] <- median.pRep.VR.S <- median(pRep.VR.S$Sensitivity)
sum.vr.sensitivity[1,3] <- mean.pRep.VR.S <- mean(pRep.VR.S$Sensitivity)
sum.vr.sensitivity[1,4] <- sd.pRep.VR.S <- sd(pRep.VR.S$Sensitivity)

# Clutch Size
sum.vr.sensitivity[2,2] <- median.mean.CS.VR.S <- median(mean.CS.VR.S$Sensitivity)
sum.vr.sensitivity[2,3] <- mean.mean.CS.VR.S <- mean(mean.CS.VR.S$Sensitivity)
sum.vr.sensitivity[2,4] <- sd.mean.CS.VR.S <- sd(mean.CS.VR.S$Sensitivity)

# Clutch Survival
sum.vr.sensitivity[3,2] <- median.S_C.VR.S <- median(S_C.VR.S$Sensitivity)
sum.vr.sensitivity[3,3] <- mean.S_C.VR.S <- mean(S_C.VR.S$Sensitivity)
sum.vr.sensitivity[3,4] <- sd.S_C.VR.S <- sd(S_C.VR.S$Sensitivity)

# Female Chick survival
sum.vr.sensitivity[4,2] <- median.sF_NB_Ch.VR.S <- median(sF_NB_Ch.VR.S$Sensitivity)
sum.vr.sensitivity[4,3] <- mean.sF_NB_Ch.VR.S <- mean(sF_NB_Ch.VR.S$Sensitivity)
sum.vr.sensitivity[4,4] <- sd.sF_NB_Ch.VR.S <- sd(sF_NB_Ch.VR.S$Sensitivity)

# Female Yearling Non-breeding Survival
sum.vr.sensitivity[5,2] <- median.sF_NB_1yr.VR.S <- median(sF_NB_1yr.VR.S$Sensitivity)
sum.vr.sensitivity[5,3] <- mean.sF_NB_1yr.VR.S <- mean(sF_NB_1yr.VR.S$Sensitivity)
sum.vr.sensitivity[5,4] <- sd.sF_NB_1yr.VR.S <- sd(sF_NB_1yr.VR.S$Sensitivity)

# Female 2 Years Non-breeding Survival
sum.vr.sensitivity[6,2] <- median.sF_NB_2yr.VR.S <- median(sF_NB_2yr.VR.S$Sensitivity)
sum.vr.sensitivity[6,3] <- mean.sF_NB_2yr.VR.S <- mean(sF_NB_2yr.VR.S$Sensitivity)
sum.vr.sensitivity[6,4] <- sd.sF_NB_2yr.VR.S <- sd(sF_NB_2yr.VR.S$Sensitivity)

# Female 3 Years Non-breeding Survival
sum.vr.sensitivity[7,2] <- median.sF_NB_3yr.VR.S <- median(sF_NB_3yr.VR.S$Sensitivity)
sum.vr.sensitivity[7,3] <- mean.sF_NB_3yr.VR.S <- mean(sF_NB_3yr.VR.S$Sensitivity)
sum.vr.sensitivity[7,4] <- sd.sF_NB_3yr.VR.S<- sd(sF_NB_3yr.VR.S$Sensitivity)

# Female Juvenile survival
sum.vr.sensitivity[8,2] <- median.sF_BN_Ju.VR.S <- median(sF_BN_Ju.VR.S$Sensitivity)
sum.vr.sensitivity[8,3] <- mean.sF_BN_Ju.VR.S <- mean(sF_BN_Ju.VR.S$Sensitivity)
sum.vr.sensitivity[8,4] <- sd.sF_BN_Ju.VR.S <- sd(sF_BN_Ju.VR.S$Sensitivity)

# Female Yearling Breeding Survival
sum.vr.sensitivity[9,2] <- median.sF_BN_1yr.VR.S <- median(sF_BN_1yr.VR.S$Sensitivity)
sum.vr.sensitivity[9,3] <- mean.sF_BN_1yr.VR.S <- mean(sF_BN_1yr.VR.S$Sensitivity)
sum.vr.sensitivity[9,4] <- sd.sF_BN_1yr.VR.S <- sd(sF_BN_1yr.VR.S$Sensitivity)

# Female 2 Years Breeding Survival
sum.vr.sensitivity[10,2] <- median.sF_BN_2yr.VR.S <- median(sF_BN_2yr.VR.S$Sensitivity)
sum.vr.sensitivity[10,3] <- mean.sF_BN_2yr.VR.S <- mean(sF_BN_2yr.VR.S$Sensitivity)
sum.vr.sensitivity[10,4] <- sd.sF_BN_2yr.VR.S <- sd(sF_BN_2yr.VR.S$Sensitivity)

# Female 3 Years Breeding Survival
sum.vr.sensitivity[11,2] <- median.sF_BN_3yr.VR.S <- median(sF_BN_3yr.VR.S$Sensitivity)
sum.vr.sensitivity[11,3] <- mean.sF_BN_3yr.VR.S <- mean(sF_BN_3yr.VR.S$Sensitivity)
sum.vr.sensitivity[11,4] <- sd.sF_BN_3yr.VR.S <- sd(sF_BN_3yr.VR.S$Sensitivity)

sum.vr.sensitivity



# Calculate median, mean and sd of elasticity

# Breeding probability
sum.vr.elasticity[1,2] <- median.pRep.VR.E <- median(pRep.VR.E$Elasticity)
sum.vr.elasticity[1,3] <- mean.pRep.VR.E <- mean(pRep.VR.E$Elasticity)
sum.vr.elasticity[1,4] <- sd.pRep.VR.E <- sd(pRep.VR.E$Elasticity)

# Clutch Size
sum.vr.elasticity[2,2] <- median.mean.CS.VR.E <- median(mean.CS.VR.E$Elasticity)
sum.vr.elasticity[2,3] <- mean.mean.CS.VR.E <- mean(mean.CS.VR.E$Elasticity)
sum.vr.elasticity[2,4] <- sd.mean.CS.VR.E <- sd(mean.CS.VR.E$Elasticity)

# Clutch Survival 
sum.vr.elasticity[3,2] <- median.S_C.VR.E <- median(S_C.VR.E$Elasticity)
sum.vr.elasticity[3,3] <- mean.S_C.VR.E <- mean(S_C.VR.E$Elasticity)
sum.vr.elasticity[3,4] <- sd.S_C.VR.E <- sd(S_C.VR.E$Elasticity)

# Female Chick survival
sum.vr.elasticity[4,2] <- median.sF_NB_Ch.VR.E <- median(sF_NB_Ch.VR.E$Elasticity)
sum.vr.elasticity[4,3] <- mean.sF_NB_Ch.VR.E <- mean(sF_NB_Ch.VR.E$Elasticity)
sum.vr.elasticity[4,4] <- sd.sF_NB_Ch.VR.E <- sd(sF_NB_Ch.VR.E$Elasticity)

# Female Yearling Non-breeding Survival
sum.vr.elasticity[5,2] <- median.sF_NB_1yr.VR.E <- median(sF_NB_1yr.VR.E$Elasticity)
sum.vr.elasticity[5,3] <- mean.sF_NB_1yr.VR.E <- mean(sF_NB_1yr.VR.E$Elasticity)
sum.vr.elasticity[5,4] <- sd.sF_NB_1yr.VR.E <- sd(sF_NB_1yr.VR.E$Elasticity)

# Female 2 Years Non-breeding Survival
sum.vr.elasticity[6,2] <- median.sF_NB_2yr.VR.E <- median(sF_NB_2yr.VR.E$Elasticity)
sum.vr.elasticity[6,3] <- mean.sF_NB_2yr.VR.E <- mean(sF_NB_2yr.VR.E$Elasticity)
sum.vr.elasticity[6,4] <- sd.sF_NB_2yr.VR.E <- sd(sF_NB_2yr.VR.E$Elasticity)

# Female 3 Years Non-breeding Survival
sum.vr.elasticity[7,2] <- median.sF_NB_3yr.VR.E <- median(sF_NB_3yr.VR.E$Elasticity)
sum.vr.elasticity[7,3] <- mean.sF_NB_3yr.VR.E <- mean(sF_NB_3yr.VR.E$Elasticity)
sum.vr.elasticity[7,4] <- sd.sF_NB_3yr.VR.E <- sd(sF_NB_3yr.VR.E$Elasticity)

# Female Juvenile survival
sum.vr.elasticity[8,2] <- median.sF_BN_Ju.VR.E <- median(sF_BN_Ju.VR.E$Elasticity)
sum.vr.elasticity[8,3] <- mean.sF_BN_Ju.VR.E <- mean(sF_BN_Ju.VR.E$Elasticity)
sum.vr.elasticity[8,4] <- sd.sF_BN_Ju.VR.E <- sd(sF_BN_Ju.VR.E$Elasticity)

# Female Yearling Breeding Survival
sum.vr.elasticity[9,2] <- median.sF_BN_1yr.VR.E <- median(sF_BN_1yr.VR.E$Elasticity)
sum.vr.elasticity[9,3] <- mean.sF_BN_1yr.VR.E <- mean(sF_BN_1yr.VR.E$Elasticity)
sum.vr.elasticity[9,4] <- sd.sF_BN_1yr.VR.E <- sd(sF_BN_1yr.VR.E$Elasticity)

# Female 2 Years Breeding Survival
sum.vr.elasticity[10,2] <- median.sF_BN_2yr.VR.E <- median(sF_BN_2yr.VR.S$Sensitivity)
sum.vr.elasticity[10,3] <- mean.sF_BN_2yr.VR.E <- mean(sF_BN_2yr.VR.S$Sensitivity)
sum.vr.elasticity[10,4] <- sd.sF_BN_2yr.VR.E <- sd(sF_BN_2yr.VR.S$Sensitivity)

# Female 3 Years Breeding Survival
sum.vr.elasticity[11,2] <- median.sF_BN_3yr.VR.E <- median(sF_BN_3yr.VR.E$Elasticity)
sum.vr.elasticity[11,3] <- mean.sF_BN_3yr.VR.E <- mean(sF_BN_3yr.VR.E$Elasticity)
sum.vr.elasticity[11,4] <- sd.sF_BN_3yr.VR.E <- sd(sF_BN_3yr.VR.E$Elasticity)

sum.vr.elasticity


# Visualization

# Prepare matrix for store results
st.vr.sen <- st.vr.elas <- matrix(data = NA, nrow = 11, ncol = 3)
colnames(st.vr.sen) <- colnames(st.vr.elas) <-  c("Median", "Mean", "SD")
rownames(st.vr.sen) <- rownames(st.vr.elas) <- c("Breeding probability", "Clutch Size", "Clutch Survival",
                                                 "Female Chick survival", 
                                                 "Female Yearling Non-breeding Survival", 
                                                 "Female 2 Years Non-breeding Survival", 
                                                 "Female 3 Years Non-breeding Survival", 
                                                 "Female Juvenile survival",
                                                 "Female Yearling Breeding Survival", 
                                                 "Female 2 Years Breeding Survival",
                                                 "Female 3 Years Breeding Survival")


# Insert data into each matrix

st.vr.sen[,1] <- as.numeric(sum.vr.sensitivity[c(1:11), 2])
st.vr.sen[,2] <- as.numeric(sum.vr.sensitivity[c(1:11), 3])
st.vr.sen[,3] <- as.numeric(sum.vr.sensitivity[c(1:11), 4])

st.vr.elas[,1] <- as.numeric(sum.vr.elasticity[c(1:11), 2])
st.vr.elas[,2] <- as.numeric(sum.vr.elasticity[c(1:11), 3])
st.vr.elas[,3] <- as.numeric(sum.vr.elasticity[c(1:11), 4])


# Save data out put
saveRDS(st.vr.sen, file = "vr.sensitivity.table.rds")
saveRDS(st.vr.elas, file = "vr.elasticity.table.rds")


# Try heat map
require(graphics); require(grDevices)


heatmap(st.vr.sen,
        col = topo.colors(40),
        Rowv = NA, Colv = NA, scale = "column",
        RowSideColors = rainbow(nrow(st.vr.sen), start = 0, end = .6),
        ColSideColors = rainbow(ncol(st.vr.sen), start = 0, end = .6),
        margins = c(5,10),
        # xlab = "Vital rate",
        # ylab = "Summary",
        main = "Sensitivity of vital rates")


heatmap(st.vr.elas, Rowv = NA, Colv = NA, scale = "column",
        RowSideColors = rainbow(nrow(st.vr.elas), start = 0, end = .2), 
        ColSideColors = rainbow(ncol(st.vr.elas), start = 0, end = .2),
        margins = c(5,10),
        # xlab = "Vital rate",
        # ylab = "Summary",
        main = "Elasticity of vital rates")



