library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Model Comparison
# Age Classes in difference season


IPM.Baseline <- readRDS("GPIPM_ChJu5060.rds")
IPM.OP1_0.1 <- readRDS("GPIPM_30yOp1.1.rds")
IPM.OP1_0.2 <- readRDS("GPIPM_30yOp1.2.rds")
IPM.OP2_0.1 <- readRDS("GPIPM_30yOp2.1.rds")
IPM.OP2_0.2 <- readRDS("GPIPM_30yOp2.2.rds")


str(IPM.Baseline)
out.mat <- list(Base = as.matrix(IPM.Baseline),
                IPM.OP1_0.1 = as.matrix(IPM.OP1_0.1),
                IPM.OP1_0.2 = as.matrix(IPM.OP1_0.2),
                IPM.OP2_0.1 = as.matrix(IPM.OP2_0.1),
                IPM.OP2_0.2 = as.matrix(IPM.OP2_0.2)
)


## Set the range of study years

StudyYears <- 2019:2051



## Extract Tmax 

ModelTmax <- length(StudyYears)


## Make an empty list for storing estimates for each model

sum.data <- list(
  Base = NA,
  Option1.0.1 = NA,
  Option1.0.2 = NA,
  Option2.0.1 = NA,
  Option2.0.2 = NA
)


## Define a short-cut function for posterior summaries (median, 95% CI)

sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.975)))
}

## For each model, extract relevant measures

ModelIDs <- c("Baseline", "Increase all Survival 10%", "Increase all Survival 20%",
              "Increase all Survival + Reproduction 10%","Increase all Survival + Reproduction 20%"
)

sum.list_VRs <- list()

# Total Population
for(i in 1:length(ModelIDs)){
  
  
  # Set sample matrix
  
  sam.mat <- out.mat[[i]]
  
  
  # Prepare matrices to store results
  Fec <- S_C <- rho <- CS <- pRep <- 
  SFemale_1_Breed <- SFemale_2_Breed <- SFemale_3_Breed <- SFemale_4_Breed <- 
  SFemale_1_Non <- SFemale_2_Non <- SFemale_3_Non <- SFemale_4_Non <-
  SMale_1_Breed <- SMale_2_Breed <- SMale_3_Breed <- SMale_4_Breed <- 
  SMale_1_Non <- SMale_2_Non <- SMale_3_Non <- SMale_4_Non <- 
  AnnSF_1 <- AnnSF_2 <- AnnSF_3 <- AnnSF_4 <- 
  AnnSM_1 <- AnnSM_2 <- AnnSM_3 <- AnnSM_4 <- 
  matrix(NA, nrow = ModelTmax, ncol = 3, dimnames = list(NULL, c('lCI', 'Median', 'uCI')))
  
  
  
  for(t in 1:ModelTmax){
    
    # Calclulate Number of each age classes 
    
    Clutchsize <- sam.mat[,paste0('CS[', t, ']')]
    Fecundity <- sam.mat[,paste0('Fec[', t, ']')]
    Broodsize <- sam.mat[,paste0('rho[', t, ']')]
    Clurchsur<- sam.mat[,paste0('S_C[', t, ']')]
    Breedprob<- sam.mat[,paste0('pRep[', t, ']')]
    SF1_Breed <- sam.mat[,paste0('sF_BN[', 1, ", ", t, ']')]
    SF2_Breed <- sam.mat[,paste0('sF_BN[', 2, ", ", t, ']')]
    SF3_Breed <- sam.mat[,paste0('sF_BN[', 3, ", ", t, ']')]
    SF4_Breed <- sam.mat[,paste0('sF_BN[', 4, ", ", t, ']')]
    SF1_Non <- sam.mat[,paste0('sF_NB[', 1, ", ", t, ']')]
    SF2_Non <- sam.mat[,paste0('sF_NB[', 2, ", ", t, ']')]
    SF3_Non <- sam.mat[,paste0('sF_NB[', 3, ", ", t, ']')]
    SF4_Non <- sam.mat[,paste0('sF_NB[', 4, ", ", t, ']')]
    SM1_Breed <- sam.mat[,paste0('sM_BN[', 1, ", ", t, ']')]
    SM2_Breed <- sam.mat[,paste0('sM_BN[', 2, ", ", t, ']')]
    SM3_Breed <- sam.mat[,paste0('sM_BN[', 3, ", ", t, ']')]
    SM4_Breed <- sam.mat[,paste0('sM_BN[', 4, ", ", t, ']')]
    SM1_Non <- sam.mat[,paste0('sM_NB[', 1, ", ", t, ']')]
    SM2_Non <- sam.mat[,paste0('sM_NB[', 2, ", ", t, ']')]
    SM3_Non <- sam.mat[,paste0('sM_NB[', 3, ", ", t, ']')]
    SM4_Non <- sam.mat[,paste0('sM_NB[', 4, ", ", t, ']')]
    
    CS[t,] <- sam.summary(Clutchsize)
    Fec[t,] <- sam.summary(Fecundity)
    rho[t,] <- sam.summary(Broodsize)
    S_C[t,] <- sam.summary(Clurchsur)
    pRep[t,] <- sam.summary(Breedprob)
      
    SFemale_1_Breed[t,] <- sam.summary(SF1_Breed)
    SFemale_2_Breed[t,] <- sam.summary(SF2_Breed)
    SFemale_3_Breed[t,] <- sam.summary(SF3_Breed)
    SFemale_4_Breed[t,] <- sam.summary(SF4_Breed)
    SFemale_1_Non[t,] <- sam.summary(SF1_Non)
    SFemale_2_Non[t,] <- sam.summary(SF2_Non)
    SFemale_3_Non[t,] <- sam.summary(SF3_Non)
    SFemale_4_Non[t,] <- sam.summary(SF4_Non)
    SMale_1_Breed[t,] <- sam.summary(SM1_Breed)
    SMale_2_Breed[t,] <- sam.summary(SM2_Breed)
    SMale_3_Breed[t,] <- sam.summary(SM3_Breed)
    SMale_4_Breed[t,] <- sam.summary(SM4_Breed)
    SMale_1_Non[t,] <- sam.summary(SM1_Non)
    SMale_2_Non[t,] <- sam.summary(SM2_Non)
    SMale_3_Non[t,] <- sam.summary(SM3_Non)
    SMale_4_Non[t,] <- sam.summary(SM4_Non)
    
    AnnSF_1[t,] <- sam.summary(SF1_Breed * SF1_Non)
    AnnSF_2[t,] <- sam.summary(SF2_Breed * SF2_Non)
    AnnSF_3[t,] <- sam.summary(SF3_Breed * SF3_Non)
    AnnSF_4[t,] <- sam.summary(SF4_Breed * SF4_Non)
    AnnSM_1[t,] <- sam.summary(SM1_Breed * SM1_Non)
    AnnSM_2[t,] <- sam.summary(SM2_Breed * SM2_Non)
    AnnSM_3[t,] <- sam.summary(SM3_Breed * SM3_Non)
    AnnSM_4[t,] <- sam.summary(SM4_Breed * SM4_Non)
    
  }
  
  model.sum.data <- data.frame(
    ModelID = ModelIDs[i],
    Year = StudyYears,
    VitalRate = rep(c('Fecundity','ClutchSurvival','BroodSize','ClutchSize',
                      'BreedingProbability','FS1BN','FS1NB','FS2BN','FS2NB',
                      'FS3BN','FS3NB','FS4BN','FS4NB','MS1BN','MS1NB','MS2BN',
                      'MS2NB','MS3BN','MS3NB','MS4BN','MS4NB','AnnFS1',
                      'AnnFS2','AnnFS3','AnnFS4','AnnMS1','AnnMS2','AnnMS3',
                      'AnnMS4'), each = ModelTmax)
    # Season = rep(c('Breeding','Non-breeding','Breeding','Non-breeding',
    #                'Breeding','Non-breeding','Breeding','Non-breeding',
    #                'Breeding','Non-breeding','Breeding','Non-breeding',
    #                'Breeding','Non-breeding','Breeding','Non-breeding',
    #                'Breeding','Non-breeding'),each = ModelTmax)
  )
  
  
  model.sum.data <- cbind(
    model.sum.data,
    rbind(Fec, S_C, rho, CS, pRep,
          SFemale_1_Breed, SFemale_1_Non, SFemale_2_Breed, SFemale_2_Non, 
          SFemale_3_Breed, SFemale_3_Non, SFemale_4_Breed, SFemale_4_Non,
          SMale_1_Breed, SMale_1_Non, SMale_2_Breed, SMale_2_Non,
          SMale_3_Breed, SMale_3_Non, SMale_4_Breed, SMale_4_Non,
          AnnSF_1, AnnSF_2, AnnSF_3, AnnSF_4, 
          AnnSM_1, AnnSM_2, AnnSM_3, AnnSM_4))
  
  
  # Store summary data
  sum.list_VRs[[i]] <- model.sum.data
  
} 

sum.list_VRs


## Make a combined data frame with all populations
allModel.data_VRs <- dplyr::bind_rows(sum.list_VRs, .id = "column_label")



# PLOT VITALRATE #
#------------#


# Take subset of data (relevant years, total numbers)
allModel.data_VRs <- subset(allModel.data_VRs, Year %in% c(2019:2032) & 
                            VitalRate %in% c('Fecundity','ClutchSurvival','BroodSize',
                                             'ClutchSize','BreedingProbability',
                                             'AnnFS1','AnnFS2','AnnFS3','AnnFS4',
                                             'AnnMS1','AnnMS2','AnnMS3','AnnMS4'))


# Add additional labels for separating scenarios
allModel.data_VRs$Effect <- dplyr::case_when(allModel.data_VRs$ModelID == "Baseline" ~ "Baseline",
                                     allModel.data_VRs$ModelID %in% c("Increase all Survival 10%", "Increase all Survival 20%") ~ "Survival",
                                     TRUE ~ "Survival & Reproduction")

allModel.data_VRs$Increase <- dplyr::case_when(allModel.data_VRs$ModelID == "Baseline" ~ "None",
                                       allModel.data_VRs$ModelID %in% c("Increase all Survival 10%", "Increase all Survival + Reproduction 10%") ~ "10%",
                                       TRUE ~ "20%")
allModel.data_VRs$Increase <- factor(allModel.data_VRs$Increase, levels = c("None", "10%", "20%"))


# Plot scenarios - Natural scale
ggplot(allModel.data_VRs, aes(x = Year, y = Median)) + 
  geom_line(aes(color = Increase, linetype = Effect)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Increase, linetype = Effect), alpha = 0.1) + 
  geom_vline(xintercept = 2026, linetype = "dashed",
             color = "black", size = 0.5) +  # After Pertubation
  geom_vline(xintercept = 2022, linetype = "dotted",
             color = "black", size = 0.5) +  # After data collecting
  # geom_text(aes(x=2020, label = "Data Collection Period"),
  #           colour="black", size = 2.5, alpha = 0.5) +
  # geom_text(aes(x=2028, label = "Management Period"),
  #           colour="black", size = 2.5,  alpha = 0.5) +
  # geom_rect(aes(xmin = 2019, xmax = 2022, fill = "green"),
  #           ymin = -Inf, ymax = Inf, alpha = 0.9,
  #           # data = presidential
  # ) +
  # scale_x_continuous(breaks = scales::breaks_width(1)) +
  ggtitle('Posterior vital rates') + 
  facet_wrap(~ VitalRate, ncol = 4, scales = "free") +
  # scale_x_continuous(breaks = scales::breaks_width(1)) +
  # scale_y_continuous(breaks = scales::breaks_width(100)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                        plot.title = element_text(face = 'bold'))


# Plot scenarios - Log scale
ggplot(allModel.data_VRs, aes(x = Year, y = log(Median))) + 
  geom_line(aes(color = Increase, linetype = Effect)) + 
  geom_ribbon(aes(ymin = log(lCI), ymax = log(uCI), fill = Increase, linetype = Effect), alpha = 0.1) + 
  # geom_vline(xintercept = 2026, linetype = "dashed",
  #            color = "black", size = 0.5) +  # After Pertubation
  # geom_vline(xintercept = 2022, linetype = "dotted",
  #            color = "black", size = 0.5) +  # After data collecting
  # geom_text(aes(x=2020, label = "Data Collection Period", y= 450),
  #           colour="black", size = 2.5, alpha = 0.5) +
  # geom_text(aes(x=2028, label = "Management Period", y= 450),
  #           colour="black", size = 2.5,  alpha = 0.5) +
  # scale_x_continuous(breaks = scales::breaks_width(1)) +
  facet_wrap(~ VitalRate, ncol = 4, scales = "free") +
  # scale_x_continuous(breaks = scales::breaks_width(1)) +
  # scale_y_continuous(breaks = scales::breaks_width(100)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                          plot.title = element_text(face = 'bold'))









