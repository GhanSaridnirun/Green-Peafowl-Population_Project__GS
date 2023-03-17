########################
# GREEN PEAFOWL IPM    #
# SENSITIVITY ANALYSIS #
########################


# Read data from rds
IPM.rhoDeriv <- readRDS("GPIPM_TwoSex_Matrix_Clutch_BreedProb_rhoDeriv.rds")

# Set as matrix
mat.ipm <- as.matrix(IPM.rhoDeriv)


## Set vital rate parameters
VR.params <- c("mean.rho", "mean.S_C", "mean.CS", 
               paste0("sF_NB[", 1:4, "]"), paste0("sF_BN[", 1:4, "]"),
               paste0("sM_NB[", 1:4, "]"), paste0("sM_BN[", 1:4, "]"))

## Re-organize data for whole posteriors of vital rate parameters
post.data <- data.frame()
out.data <- melt(mat.ipm[,VR.params])
colnames(out.data) <- c('Sample', 'Parameter', 'Estimate')
post.data <- rbind(post.data, out.data)









# 1. Set parameter values #
#-------------------------#



# 2. Build projection matrix #
#----------------------------#



# 3. Calculate asymptotic population growth rate #
#------------------------------------------------#



# 4. Calculate matrix element sensitivities #
#-------------------------------------------#



# 5. Calculate matrix element elasticities #
#------------------------------------------#



# 6. Calculate vital rate sensitivities #
#---------------------------------------#



# 7. Calculate vital rate elasticieites #
#---------------------------------------#

