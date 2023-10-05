library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Model Comparison
# Age Classes in difference season

IPM.CJ4050 <- readRDS("GPIPM_ChJu4050.rds")
IPM.CJ5060 <- readRDS("GPIPM_ChJu5060.rds")
IPM.CJ6070 <- readRDS("GPIPM_ChJu6070.rds")


out.mat <- list(IPM.CJ4050 = as.matrix(IPM.CJ4050),
                IPM.CJ5060 = as.matrix(IPM.CJ5060),
                IPM.CJ6070 = as.matrix(IPM.CJ6070)
)


## Set the range of study years

StudyYears <- 2019:2052



## Extract Tmax 

ModelTmax <- length(StudyYears)


## Make an empty list for storing estimates for each model

sum.data <- list(
  IPM.CJ4050 = NA,
  IPM.CJ5060 = NA,
  IPM.CJ6070 = NA
)


## Define a short-cut function for posterior summaries (median, 95% CI)

sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.975)))
}



## For each model, extract relevant measures

ModelIDs <- c("Low(40-50)", "ModelSetUp(50-60)", "High(60-70)"
)

sum.list_BN_NB <- list()

# Total Population
for(i in 1:length(ModelIDs)){
  
  
  # Set sample matrix
  
  sam.mat <- out.mat[[i]]
  
  
  # Prepare matrices to store results
  NFemale_1_Breed <- NFemale_2_Breed <- NFemale_3_Breed <- NFemale_4_Breed <- 
    NFemale_1_Non <- NFemale_2_Non <- NFemale_3_Non <- NFemale_4_Non <-
    NMale_1_Breed <- NMale_2_Breed <- NMale_3_Breed <- NMale_4_Breed <- 
    NMale_1_Non <- NMale_2_Non <- NMale_3_Non <- NMale_4_Non <- 
    NTot_Breed <- NTot_Non <- 
    matrix(NA, nrow = ModelTmax, ncol = 3, dimnames = list(NULL, c('lCI', 'Median', 'uCI')))
  
  
  
  for(t in 1:ModelTmax){
    
    # Calclulate Number of each age classes 
    
    NF1_Breed <- sam.mat[,paste0('NBreedF[', 1, ", ", t, ']')]
    NF2_Breed <- sam.mat[,paste0('NBreedF[', 2, ", ", t, ']')]
    NF3_Breed <- sam.mat[,paste0('NBreedF[', 3, ", ", t, ']')]
    NF4_Breed <- sam.mat[,paste0('NBreedF[', 4, ", ", t, ']')]
    NF1_Non <- sam.mat[,paste0('NNonF[', 1, ", ", t, ']')]
    NF2_Non <- sam.mat[,paste0('NNonF[', 2, ", ", t, ']')]
    NF3_Non <- sam.mat[,paste0('NNonF[', 3, ", ", t, ']')]
    NF4_Non <- sam.mat[,paste0('NNonF[', 4, ", ", t, ']')]
    NM1_Breed <- sam.mat[,paste0('NBreedM[', 1, ", ", t, ']')]
    NM2_Breed <- sam.mat[,paste0('NBreedM[', 2, ", ", t, ']')]
    NM3_Breed <- sam.mat[,paste0('NBreedM[', 3, ", ", t, ']')]
    NM4_Breed <- sam.mat[,paste0('NBreedM[', 4, ", ", t, ']')]
    NM1_Non <- sam.mat[,paste0('NNonM[', 1, ", ", t, ']')]
    NM2_Non <- sam.mat[,paste0('NNonM[', 2, ", ", t, ']')]
    NM3_Non <- sam.mat[,paste0('NNonM[', 3, ", ", t, ']')]
    NM4_Non <- sam.mat[,paste0('NNonM[', 4, ", ", t, ']')]
    
    
    
    NFemale_1_Breed[t,] <- sam.summary(NF1_Breed)
    NFemale_2_Breed[t,] <- sam.summary(NF2_Breed)
    NFemale_3_Breed[t,] <- sam.summary(NF3_Breed)
    NFemale_4_Breed[t,] <- sam.summary(NF4_Breed)
    NFemale_1_Non[t,] <- sam.summary(NF1_Non)
    NFemale_2_Non[t,] <- sam.summary(NF2_Non)
    NFemale_3_Non[t,] <- sam.summary(NF3_Non)
    NFemale_4_Non[t,] <- sam.summary(NF4_Non)
    NMale_1_Breed[t,] <- sam.summary(NM1_Breed)
    NMale_2_Breed[t,] <- sam.summary(NM2_Breed)
    NMale_3_Breed[t,] <- sam.summary(NM3_Breed)
    NMale_4_Breed[t,] <- sam.summary(NM4_Breed)
    NMale_1_Non[t,] <- sam.summary(NM1_Non)
    NMale_2_Non[t,] <- sam.summary(NM2_Non)
    NMale_3_Non[t,] <- sam.summary(NM3_Non)
    NMale_4_Non[t,] <- sam.summary(NM4_Non)
    
    NTot_Breed[t,] <- sam.summary(NF1_Breed + NF2_Breed + NF3_Breed + NF4_Breed + 
                                    NM1_Breed + NM2_Breed + NM3_Breed + NM4_Breed)
    
    NTot_Non[t,] <- sam.summary(NF1_Non + NF2_Non + NF3_Non + NF4_Non + 
                                  NM1_Non + NM2_Non + NM3_Non + NM4_Non)
    
    
  }
  
  model.sum.data <- data.frame(
    ModelID = ModelIDs[i],
    Year = StudyYears,
    AgeClass = rep(c('Female1_BN', 'Female1_NB','Female2_BN', 'Female2_NB',
                     'Female3_BN', 'Female3_NB','Female4_BN', 'Female4_NB',
                     'Male1_BN', 'Male1_NB','Male2_BN', 'Male2_NB',
                     'Male3_BN', 'Male3_NB','Male4_BN', 'Male4_NB',
                     'Total_BN', 'Total_NB'), each = ModelTmax),
    Season = rep(c('Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding', 'Non-breeding'), each = ModelTmax)
  )
  
  
  model.sum.data <- cbind(
    model.sum.data,
    rbind(NFemale_1_Breed, NFemale_1_Non, NFemale_2_Breed, NFemale_2_Non, 
          NFemale_3_Breed, NFemale_3_Non, NFemale_4_Breed, NFemale_4_Non,
          NMale_1_Breed, NMale_1_Non, NMale_2_Breed, NMale_2_Non,
          NMale_3_Breed, NMale_3_Non, NMale_4_Breed, NMale_4_Non,
          NTot_Breed, NTot_Non))
  
  
  
  # Store summary data
  sum.list_BN_NB[[i]] <- model.sum.data
  
} 

sum.list_BN_NB


## Make a combined data frame with all populations
allModel.data_BN_NB <- dplyr::bind_rows(sum.list_BN_NB, .id = "column_label")

# PLOT TOTAL #
#------------#

# Take subset of data (relevant years, total numbers)
TotalData <- subset(allModel.data_BN_NB, Year %in% c(2022:2035) & AgeClass %in% c('Total_BN', 'Total_NB'))

# # Drop first year non-breeding season (not estimated)
# TotalData <- subset(TotalData, !(Season == "Non-breeding" & Year == 2019))
# 
# # Add additional labels for separating scenarios
# TotalData$Effect <- dplyr::case_when(TotalData$ModelID == "Baseline" ~ "Baseline",
#                                      TotalData$ModelID %in% c("Increase all Survival 10%", "Increase all Survival 20%") ~ "Survival",
#                                      TRUE ~ "Survival & Reproduction")
# 
# TotalData$Increase <- dplyr::case_when(TotalData$ModelID == "Baseline" ~ "None",
#                                        TotalData$ModelID %in% c("Increase all Survival 10%", "Increase all Survival + Reproduction 10%") ~ "10%",
#                                        TRUE ~ "20%")
TotalData$Increase <- factor(TotalData$ModelID, levels = c("ModelSetUp(50-60)", 
                                                           "Low(40-50)",
                                                           "High(60-70)"))

# Plot scenarios - Natural scale
ggplot(TotalData, aes(x = Year, y = Median)) + 
  geom_line(aes(color = ModelID, linetype = ModelID)) + 
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = ModelID, linetype = ModelID), alpha = 0.1) + 
  # geom_vline(xintercept = 2026, linetype = "dashed",
  #            color = "black", size = 0.5) +  # After Pertubation
  # geom_vline(xintercept = 2022, linetype = "dotted",
  #            color = "black", size = 0.5) +  # After data collecting
  # geom_text(aes(x=2020, label = "Data Collection Period", y= 450),
  #           colour="black", size = 2.5, alpha = 0.5) +
  # geom_text(aes(x=2028, label = "Management Period", y= 450),
  #           colour="black", size = 2.5,  alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "black", size = 1) +  
  facet_wrap(~ Season, ncol = 1) +
  ylab("Total Population Size") + 
  ggtitle('Comparison of predicted population size from difference chick and juvenile survival rates') + 
  scale_color_manual(values = magma(4)[1:3]) + 
  scale_fill_manual(values = magma(4)[1:3]) + 
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                          plot.title = element_text(face = 'bold'))

# # Plot scenarios - Log scale
# ggplot(TotalData, aes(x = Year, y = Median)) + 
#   geom_line(aes(color = ModelID, linetype = ModelID)) + 
#   geom_ribbon(aes(ymin = log(lCI), ymax = log(uCI), fill = ModelID, linetype = ModelID), alpha = 0.1) + 
#   # geom_vline(xintercept = 2026, linetype = "dashed",
#   #            color = "black", size = 0.5) +  # After Pertubation
#   # geom_vline(xintercept = 2022, linetype = "dotted",
#   #            color = "black", size = 0.5) +  # After data collecting
#   # geom_text(aes(x=2020, label = "Data Collection Period", y= 450),
#   #           colour="black", size = 2.5, alpha = 0.5) +
#   # geom_text(aes(x=2028, label = "Management Period", y= 450),
#   #           colour="black", size = 2.5,  alpha = 0.5) +
#   geom_hline(yintercept = 0, linetype = "solid",
#              color = "black", size = 1) +  
#   facet_wrap(~ Season, ncol = 1) +
#   ylab("Total Population Size") + 
#   scale_color_manual(values = magma(4)[1:3]) + 
#   scale_fill_manual(values = magma(4)[1:3]) + 
#   scale_x_continuous(breaks = scales::breaks_width(1)) +
#   scale_y_continuous(breaks = scales::breaks_width(100)) +
#   theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
#                           plot.title = element_text(face = 'bold'))


