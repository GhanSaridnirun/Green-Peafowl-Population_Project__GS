########################
# GREEN PEAFOWL IPM    #
# SENSITIVITY ANALYSIS #
########################

library(coda)

# Load Function for sensitivity and elasticity calculation
source("R/GPeafowlMatrix_SensElasAnalysis.R")

# Read data from rds
IPM.rhoDeriv <- readRDS("GPIPM_TwoSex_Matrix_Clutch_BreedProb_rhoDeriv.rds")


# Set as matrix
mat.ipm <- as.matrix(IPM.rhoDeriv)

# MCMCvis::MCMCsummary(mat.ipm)
# MCMCvis::MCMCpstr(mat.ipm)



# 1. Set parameter values #
#-------------------------#

## Productivity 

pRep <- 0.5   
mean.CS <- mat.ipm[1, "mean.CS"]
gamma <- 0.5 
S_C <- mat.ipm[1, "S_C[1]"]


## Survival

sF_NB <- c(mat.ipm[1, "sF_NB[1]"], mat.ipm[1, "sF_NB[2]"], 
           mat.ipm[1, "sF_NB[3]"], mat.ipm[1, "sF_NB[4]"])

sF_BN <- c(mat.ipm[1, "sF_BN[1]"], mat.ipm[1, "sF_BN[2]"], 
           mat.ipm[1, "sF_BN[3]"], mat.ipm[1, "sF_BN[4]"])

sM_NB <- c(mat.ipm[1, "sM_NB[1]"], mat.ipm[1, "sM_NB[2]"], 
           mat.ipm[1, "sM_NB[3]"], mat.ipm[1, "sM_NB[4]"])

sM_BN <- c(mat.ipm[1, "sM_BN[1]"], mat.ipm[1, "sM_BN[2]"], 
           mat.ipm[1, "sM_BN[3]"], mat.ipm[1, "sM_BN[4]"])



# 2. Build projection matrix #
#----------------------------#

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

## Extract the population number of each age classes

ChF_NB <- mat.ipm[1, "NNonF[1, 1]"] # Female Chick Non-breeding 
F1y_NB <- mat.ipm[1, "NNonF[2, 1]"] # Female 1 year Non-breeding 
F2y_NB <- mat.ipm[1, "NNonF[3, 1]"] # Female 2 year Non-breeding
F3y_NB <- mat.ipm[1, "NNonF[3, 1]"] # Female 3 year Non-breeding
ChM_NB <- mat.ipm[1, "NNonM[1, 1]"] # Male Chick Non-breeding 
M1y_NB <- mat.ipm[1, "NNonM[2, 1]"] # Male 1 year Non-breeding 
M2y_NB <- mat.ipm[1, "NNonM[3, 1]"] # Male 2 year Non-breeding 
M3y_NB <- mat.ipm[1, "NNonM[3, 1]"] # Male 3 year Non-breeding 

JuF_BN <- mat.ipm[1, "NBreedF[1, 1]"] # Female Juvenile Breeding 
F1y_BN <- mat.ipm[1, "NBreedF[2, 1]"] # Female 1 year Breeding 
F2y_BN <- mat.ipm[1, "NBreedF[3, 1]"] # Female 2 year Breeding 
F3y_BN <- mat.ipm[1, "NBreedF[3, 1]"] # Female 3 year Breeding 
JuM_BN <- mat.ipm[1, "NBreedM[1, 1]"] # Male Juvenile Breeding 
M1y_BN <- mat.ipm[1, "NBreedM[2, 1]"] # Male 1 year Breeding 
M2y_BN <- mat.ipm[1, "NBreedM[3, 1]"] # Male 2 year Breeding 
M3y_BN <- mat.ipm[1, "NBreedM[3, 1]"] # Male 3 year Breeding 


## Set initial population to use in simulations

N1.ann <- c(JuF_BN, JuM_BN, # Chicks
            F1y_BN, M1y_BN, # 1 yr olds
            F2y_BN, M2y_BN, # 2yr olds
            F3y_BN, M3y_BN) # 3+yr olds

N1.season <- c(N1.ann[1:2],
               0, 0,
               N1.ann[3:4],
               0, 0,
               N1.ann[5:6],
               0, 0,
               N1.ann[7:8],
               0, 0)

## Annual matrix
# Using eigen decomposition
lambda.ann <- eigen(mat.ann$A)$values[1]
lambda.ann

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



# 4. Calculate matrix element sensitivities #
#-------------------------------------------#

# Set perturbation factor
dy <- 1e-5

# Set up original and sensitivity matrices
A_orig <- mat.ann$A
sens <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig), dimnames = dimnames(A_orig))

# Calculate Sensitivities
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

# Calculate Elasticities
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

# Visualize sensitivities for all vital rates 
library(ggplot2)
ggplot(VR.sens.data) + 
  geom_bar(aes(x = VitalRate, y = Sensitivity), stat = "identity") + 
  coord_flip() + 
  theme_classic()

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

# Visualize sensitivities for all vital rates 
library(ggplot2)
ggplot(VR.elas.data) + 
  geom_bar(aes(x = VitalRate, y = Elasticity), stat = "identity") + 
  coord_flip() + 
  theme_classic()




