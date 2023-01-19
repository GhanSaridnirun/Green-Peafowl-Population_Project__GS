library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


IPM.rhoEst <- readRDS("GPIPM_TwoSex_Matrix_Clutch_BreedProb_rhoEst.rds")
IPM.rhoDeriv <- readRDS("GPIPM_TwoSex_Matrix_Clutch_BreedProb_rhoDeriv.rds")


out.mat <- list(rhoEst = as.matrix(IPM.rhoEst),
                rhoDeriv = as.matrix(IPM.rhoDeriv)
)



## Set the range of study years

StudyYears <- 2019:2072



## Extract Tmax 

ModelTmax <- length(StudyYears)


## Make an empty list for storing estimates for each model

sum.data <- list(
  rhoEst = NA,
  rhoDeriv = NA
)


## Define a short-cut function for posterior summaries (median, 95% CI)

sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.975)))
}



## For each model, extract relevant measures

ModelIDs <- c("rho estimated", "rho derived from CS")

sum.list <- list()


for(i in 1:length(ModelIDs)){
  
  
  # Set sample matrix
  
  sam.mat <- out.mat[[i]]
  
  
  # Prepare matrices to store results
  NTOT_Breed <- NTOT_Non <- matrix(NA, nrow = ModelTmax, ncol = 3, dimnames = list(NULL, c('lCI', 'Median', 'uCI')))
  
 
  for(t in 1:ModelTmax){
    
    # Calclulate NTOT in both seasons
    N_Breed <- rowSums(sam.mat[,paste0('NBreedF[', 1:4, ", ", t, ']')]) + 
               rowSums(sam.mat[,paste0('NBreedM[', 1:4, ", ", t, ']')])
    
    N_Non <- rowSums(sam.mat[,paste0('NNonF[', 1:4, ", ", t, ']')]) + 
             rowSums(sam.mat[,paste0('NNonM[', 1:4, ", ", t, ']')])
    
    
    # Make posterior summary for NTOT in each year
    NTOT_Breed[t,] <- sam.summary(N_Breed)
    NTOT_Non[t,] <- sam.summary(N_Non)
  }
  
  
  # Organise summaries in a data frame
  model.sum.data <- data.frame(
    ModelID = ModelIDs[i],
    Year = StudyYears,
    Parameter = rep(c('NTOT_Breed', 'NTOT_Non'), each = ModelTmax)
  )
  
  model.sum.data <- cbind(
    model.sum.data,
    rbind(NTOT_Breed, NTOT_Non))
  
  # Store summary data
  sum.list[[i]] <- model.sum.data
  
} 

sum.list


## Make a combined data frame with all populations
allModel.data <- dplyr::bind_rows(sum.list, .id = "column_label")



## Plot Population sizes

ggplot(allModel.data, aes(x = Year, y = Median, group = ModelID)) +
  geom_line(aes(color = ModelID)) +
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = ModelID), alpha = 0.2) +
  ggtitle('Total Population Size') + 
  ylab('Estimate') + 
  facet_wrap(~Parameter, ncol = 1, scales = 'free_y') +
  theme_bw() + theme(legend.position = "right", panel.grid.major.y = element_blank()
                   , panel.grid.minor.y = element_blank(),
                   axis.text.x = element_text(angle = 45, vjust = 0.5),
                   plot.title = element_text(face = 'bold'))




#-----------------------------------------------------------------------------


## Set vital rate parameters to plot
VR.params <- c("mean.rho", "mean.S_C", "mean.CS", 
               paste0("sF_NB[", 1:4, "]"), paste0("sF_BN[", 1:4, "]"),
               paste0("sM_NB[", 1:4, "]"), paste0("sM_BN[", 1:4, "]"))


## Re-organize data for whole posteriors of vital rate parameters
post.data <- data.frame()

for(i in 1:2){
  out.data <- melt(out.mat[[i]][,VR.params])
  colnames(out.data) <- c('Sample', 'Parameter', 'Estimate')
  out.data$ModelID <- ModelIDs[i]
  post.data <- rbind(post.data, out.data)
}



# Plotting: Vital rate posteriors, across models 

pdf('ModelComp_VRs.pdf', width = 10, height = 8)
ggplot(post.data, aes(x = Estimate, group = ModelID)) + 
  geom_density(aes(color = ModelID, fill = ModelID), alpha = 0.5) +
  facet_wrap(~Parameter, scales = 'free') +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(face = 'bold'))
dev.off()




