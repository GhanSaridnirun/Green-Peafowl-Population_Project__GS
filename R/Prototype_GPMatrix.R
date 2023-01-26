
###############
# FULL MATRIX #
###############

# Number Assignment 






N <- matrix(0, nrow = 16, ncol = 16)
dimnames(N) <- list(c("ChF_NB", "ChM_NB", "JuF_BN", "JuM_BN",
                      "1yF_NB", "1yM_NB", "1yF_BN", "1yM_BN",
                      "2yF_NB", "2yM_NB", "2yF_BN", "2yM_BN",
                      "3yF_NB", "3yM_NB", "3yF_BN", "3yM_BN"),
                    
                    c("ChF_NB", "ChM_NB", "JuF_BN", "JuM_BN",
                      "1yF_NB", "1yM_NB", "1yF_BN", "1yM_BN",
                      "2yF_NB", "2yM_NB", "2yF_BN", "2yM_BN",
                      "3yF_NB", "3yM_NB", "3yF_BN", "3yM_BN"))


pop.mat <- N
pop.mat

pop.mat.num <- pop.mat


# Row index = gets produced
# Column index = produces

# Female chicks (NB)
pop.mat["ChF_NB", "2yF_BN"] <- "pRep * mean.CS * S_C * gamma"
pop.mat["ChF_NB", "3yF_BN"] <- "pRep * mean.CS * S_C * gamma"

# Male chicks (NB)
pop.mat["ChF_NB", "2yF_BN"] <- "pRep * mean.CS * S_C * (1-gamma)"
pop.mat["ChF_NB", "3yF_BN"] <- "pRep * mean.CS * S_C * (1-gamma)"

# Female juvenile (BN)
pop.mat["JuF_BN", "ChF_NB"] <- "sF_NB[1]"

# Male juvenile (BN)
pop.mat["JuM_BN", "ChM_NB"] <- "sM_NB[1]"

# Female yearling (NB)
pop.mat["1yF_NB", "JuF_BN"] <- "sF_BN[1]"

# Male yearling (NB)
pop.mat["1yM_NB", "JuM_BN"] <- "sM_BN[1]"

# Female yearling (BN)
pop.mat["1yF_BN", "1yF_NB"] <- "sF_NB[2]"

# Male yearling (BN)
pop.mat["1yM_BN", "1yM_NB"] <- "sM_NB[2]"

# Female 2-yr (NB)
pop.mat["2yF_NB", "1yF_BN"] <- "sF_BN[2]"

# Male 2-yr (NB)
pop.mat["2yM_NB", "1yM_BN"] <- "sM_BN[2]"

# Female 2-yr (BN)
pop.mat["2yF_BN", "2yF_NB"] <- "sF_NB[3]"

# Male 2-yr (BN)
pop.mat["2yM_BN", "2yM_NB"] <- "sM_NB[3]"

# Female 3-yr (NB)
pop.mat["3yF_NB", "2yF_BN"] <- "sF_BN[3]"
pop.mat["3yF_NB", "3yF_BN"] <- "sF_BN[4]"

# Male 3-yr (NB)
pop.mat["3yM_NB", "2yM_BN"] <- "sM_BN[3]"
pop.mat["3yM_NB", "3yM_BN"] <- "sM_BN[4]"

# Females 3-yr (BN)
pop.mat["3yF_BN", "3yF_NB"] <- "sF_NB[4]"

# Males 3-yr (BN)
pop.mat["3yM_BN", "3yM_NB"] <- "sM_NB[4]"



#################
# ANNUAL MATRIX #
#################

# Census at start of non-breeding season (NB), chicks just produced

N <- matrix(0, nrow = 8, ncol = 8)
dimnames(N) <- list(c("ChF", "ChM",
                      "1yF", "1yM",
                      "2yF", "2yM",
                      "3yF", "3yM"),
                    
                    c("ChF", "ChM",
                      "1yF", "1yM",
                      "2yF", "2yM",
                      "3yF", "3yM"))


ann.pop.mat <- N
ann.pop.mat

ann.pop.mat.num <- ann.pop.mat


# Row index = gets produced
# Column index = produces

# Female chicks (NB)
ann.pop.mat["ChF", "1yF"] <- "sF_NB[2] * pRep * mean.CS * S_C * gamma"
ann.pop.mat["ChF", "2yF"] <- "sF_NB[3] * pRep * mean.CS * S_C * gamma"
ann.pop.mat["ChF", "3yF"] <- "sF_NB[4] * pRep * mean.CS * S_C * gamma"


# Male chicks (NB)
ann.pop.mat["ChM", "1yF"] <- "sF_NB[2] * pRep * mean.CS * S_C * (1-gamma)"
ann.pop.mat["ChM", "2yF"] <- "sF_NB[3] * pRep * mean.CS * S_C * (1-gamma)"
ann.pop.mat["ChM", "3yF"] <- "sF_NB[4] * pRep * mean.CS * S_C * (1-gamma)"


# Female yearling (NB)
ann.pop.mat["1yF", "ChF"] <- "sF_NB[1] * sF_BN[1]"

# Male yearling (NB)
ann.pop.mat["1yM", "ChM"] <- "sM_NB[1] * sM_BN[1]"

# Female 2-yr (NB)
ann.pop.mat["2yF", "1yF"] <- "sF_NB[2] * sF_BN[2]"

# Male 2-yr (NB)
ann.pop.mat["2yM", "1yM"] <- "sM_NB[2] * sM_BN[2]"

# Female 3-yr (NB)
ann.pop.mat["3yF", "2yF"] <- "sF_NB[3] * sF_BN[3]"
ann.pop.mat["3yF", "3yF"] <- "sF_NB[4] * sF_BN[4]"

# Male 3-yr (NB)
ann.pop.mat["3yM", "2yM"] <- "sM_NB[3] * sM_BN[3]"
ann.pop.mat["3yM", "3yM"] <- "sM_NB[4] * sM_BN[4]"


# Next tasks: 
# - Make numerical matrices
# - Calculate asymptotic population growth rate = dominant right eigenvalue of the matrix using eigen()
# - Compare pop. growth rate: should be the same


###########################
## Applying the numberic ##
###########################

# Vital rates

# Productivity 

pRep <- 0.5
mean.CS <- 5
gamma <- 0.5
S_C <- 0.5

# Survival rate 

sFCh_NB <- 0.5  # Female Chick survival Non-breeding
sMCh_NB <- 0.5  # Male Chick survival Non-breeding
sFCh_BN <- 0.5  # Female Juvenile survival Breeding
sMCh_BN <- 0.5  # Male Juvenile survival Breeding

sF1y_NB <- 0.6  # Female Yearling survival Non-breeding 
sM1y_NB <- 0.6  # Male Yearling survival Non-breeding 
sF1y_BN <- 0.6  # Female Yearling survival Breeding
sM1y_BN <- 0.6  # Male Yearling survival Breeding

sF2y_NB <- 0.7  # Female 2 years old survival Non-breeding 
sM2y_NB <- 0.7  # Male 2 years old survival Non-breeding 
sF2y_BN <- 0.7  # Female 2 years old survival Breeding
sM2y_BN <- 0.7  # Male 2 years old survival Breeding

sF3y_NB <- 0.8  # Female 3 years old survival Non-breeding 
sM3y_NB <- 0.8  # Male 3 years old survival Non-breeding 
sF3y_BN <- 0.8  # Female 3 years old survival Breeding
sM3y_BN <- 0.8  # Male 3 years old survival Breeding



N <- matrix(0, nrow = 8, ncol = 8)
dimnames(N) <- list(c("ChF", "ChM",
                      "1yF", "1yM",
                      "2yF", "2yM",
                      "3yF", "3yM"),
                    
                    c("ChF", "ChM",
                      "1yF", "1yM",
                      "2yF", "2yM",
                      "3yF", "3yM"))


ann.pop.mat <- N
ann.pop.mat

ann.pop.mat.num <- ann.pop.mat


# Row index = gets produced
# Column index = produces

# Female chicks (NB)
ann.pop.mat["ChF", "1yF"] <- sF1y_NB * pRep * mean.CS * S_C * gamma
ann.pop.mat["ChF", "2yF"] <- sF2y_NB * pRep * mean.CS * S_C * gamma
ann.pop.mat["ChF", "3yF"] <- sF3y_NB * pRep * mean.CS * S_C * gamma


# Male chicks (NB)
ann.pop.mat["ChM", "1yF"] <- sF1y_NB * pRep * mean.CS * S_C * (1-gamma)
ann.pop.mat["ChM", "2yF"] <- sF2y_NB * pRep * mean.CS * S_C * (1-gamma)
ann.pop.mat["ChM", "3yF"] <- sF3y_NB * pRep * mean.CS * S_C * (1-gamma)


# Female yearling (NB)
ann.pop.mat["1yF", "ChF"] <- sFCh_NB * sFCh_BN

# Male yearling (NB)
ann.pop.mat["1yM", "ChM"] <- sMCh_NB * sMCh_BN

# Female 2-yr (NB)
ann.pop.mat["2yF", "1yF"] <- sF1y_NB * sF1y_BN

# Male 2-yr (NB)
ann.pop.mat["2yM", "1yM"] <- sM1y_NB * sM1y_BN

# Female 3-yr (NB)
ann.pop.mat["3yF", "2yF"] <- sF2y_NB * sF2y_BN
ann.pop.mat["3yF", "3yF"] <- sF3y_NB * sF3y_BN

# Male 3-yr (NB)
ann.pop.mat["3yM", "2yM"] <- sM2y_NB * sM2y_BN
ann.pop.mat["3yM", "3yM"] <- sM3y_NB * sM3y_BN


# Calculate lambda (= dominant right eigenvalue)
eigen(ann.pop.mat)$values[1]

# Calculate stable age structure
eigen(ann.pop.mat)$vectors[,1] / sum(eigen(ann.pop.mat)$vectors[,1])

# Using popbio package
popbio::eigen.analysis(ann.pop.mat)
























