library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Model Comparison
# Age Classes in difference season

IPM.Baseline <- readRDS("GPeafowlIPM_forecastModel_30Ys_Baseline.rds")
IPM.OP1_0.1 <- readRDS("GPeafowlIPM_forecastModel_30Ys_OP01_0.1.rds")
IPM.OP1_0.2 <- readRDS("GPeafowlIPM_forecastModel_30Ys_OP01_0.2.rds")
IPM.OP2_0.1 <- readRDS("GPeafowlIPM_forecastModel_30Ys_OP02_0.1.rds")
IPM.OP2_0.2 <- readRDS("GPeafowlIPM_forecastModel_30Ys_OP02_0.2.rds")


out.mat <- list(Base = as.matrix(IPM.Baseline),
                IPM.OP1_0.1 = as.matrix(IPM.OP1_0.1),
                IPM.OP1_0.2 = as.matrix(IPM.OP1_0.2),
                IPM.OP2_0.1 = as.matrix(IPM.OP2_0.1),
                IPM.OP2_0.2 = as.matrix(IPM.OP2_0.2)
)


## Set the range of study years

StudyYears <- 2019:2049



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
    
    
    
  }
  
  model.sum.data <- data.frame(
    ModelID = ModelIDs[i],
    Year = StudyYears,
    AgeClass = rep(c('Female1_BN', 'Female1_NB','Female2_BN', 'Female2_NB',
                     'Female3_BN', 'Female3_NB','Female4_BN', 'Female4_NB',
                     'Male1_BN', 'Male1_NB','Male2_BN', 'Male2_NB',
                     'Male3_BN', 'Male3_NB','Male4_BN', 'Male4_NB'), each = ModelTmax),
    Season = rep(c('Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding','Non-breeding', 'Breeding','Non-breeding',
                   'Breeding','Non-breeding', 'Breeding','Non-breeding'),each = ModelTmax)
  )
  
  
  model.sum.data <- cbind(
    model.sum.data,
    rbind(NFemale_1_Breed, NFemale_1_Non, NFemale_2_Breed, NFemale_2_Non, 
          NFemale_3_Breed, NFemale_3_Non, NFemale_4_Breed, NFemale_4_Non,
          NMale_1_Breed, NMale_1_Non, NMale_2_Breed, NMale_2_Non,
          NMale_3_Breed, NMale_3_Non, NMale_4_Breed, NMale_4_Non))
  
  
  
  # Store summary data
  sum.list_BN_NB[[i]] <- model.sum.data
  
} 

sum.list_BN_NB


## Make a combined data frame with all populations
allModel.data_BN_NB <- dplyr::bind_rows(sum.list_BN_NB, .id = "column_label")

# Selected years period to plot
CropYear <- subset(allModel.data_BN_NB, Year %in% c(2019:2032))


CropYear$ModelID <- factor(CropYear$ModelID, levels = c("Baseline", 
                                                        "Increase all Survival 10%", 
                                                        "Increase all Survival 20%",
                                                        "Increase all Survival + Reproduction 10%",
                                                        "Increase all Survival + Reproduction 20%"))


# Plot area: Breeding
ggplot(CropYearBreeding, aes(x = Year, y = Median)) +
  geom_area(aes(fill = AgeClass), alpha = 0.9) +
  scale_fill_viridis(labels=c('Female[1]', 'Female[2]','Female[3]', 'Female[4]',
                              'Male[1]', 'Male[2]','Male[3]', 'Male[4]'), discrete = T) +
  facet_wrap(~ModelID, scales = "free_y", ncol = 1) + 
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::breaks_width(100)) +
  # guides(col = guide_legend(ncol = 2)) +
  theme_classic() + theme(legend.position = "bottom",
                          axis.text.x = element_text(angle = 45, vjust = 0.5),
                          plot.title = element_text(face = 'bold'),) +
  ggtitle('Population Size in Breeding Season') + 
  ylab('Estimate Population Size') +
  geom_vline(xintercept = 2026, linetype="dashed",
             color = "black", size=1) +  # After Pertubation
  geom_vline(xintercept = 2022, linetype="dotted",
             color = "black", size=1) +  # After data collecting
  geom_text(aes(x=2020, label="Data Collection Period", y= 340),
            colour="black", size = 2.5, angle=0, alpha = 0.5) +
  geom_text(aes(x=2028, label="Management Period", y= 340),
            colour="black", size = 2.5,  angle=0, alpha = 0.5)


# Plot area: Non-Breeding
ggplot(CropYearNonBreeding, aes(x = Year, y = Median)) +
  geom_area(aes(fill = AgeClass), alpha = 0.9) +
  scale_fill_viridis(labels=c('Female[1]', 'Female[2]','Female[3]', 'Female[4]',
                              'Male[1]', 'Male[2]','Male[3]', 'Male[4]'), discrete = T) +
  facet_wrap(~ModelID, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(breaks = scales::breaks_width(100)) + 
  # guides(col = guide_legend(ncol = 8)) +
  theme_classic() + theme(legend.position = "bottom",
                          axis.text.x = element_text(angle = 45, vjust = 0.5),
                          plot.title = element_text(face = 'bold'),) +
  ggtitle('Population Size in Non-breeding Season') + 
  ylab('Estimate Population Size') +
  geom_vline(xintercept = 2026, linetype="dashed",
             color = "black", size=1) +  # After Pertubation
  geom_vline(xintercept = 2022, linetype="dotted",
             color = "black", size=1) +  # After data collecting
  geom_text(aes(x=2020, label="Data Collection Period", y= 450),
            colour="black", size = 2.5, angle=0, alpha = 0.5) +
  geom_text(aes(x=2029, label="Management Period", y= 450),
            colour="black", size = 2.5,  angle=0, alpha = 0.5)




#-------------------------------------------------------------------------------