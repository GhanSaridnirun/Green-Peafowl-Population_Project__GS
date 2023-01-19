
# Test Run 

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)




rhopRep.IPM <- readRDS('[A__PeafowlIPM_rhoxpRep.rds')
ChSpRep.IPM <- readRDS('[A__PeafowlIPM_CSxpRep.rds')
rhoprior.IPM <- readRDS('[A__PeafowlIPM_rhoxpRep_rho_onlyprior.rds')
Choutpro.IPM <- readRDS('[A__PeafowlIPM_pRepxNoCSinProcess.rds')



out.mat <- list(rhpR = as.matrix(rhopRep.IPM),
                ChpR = as.matrix(ChSpRep.IPM),
                rhpr = as.matrix(rhoprior.IPM),
                chou = as.matrix(Choutpro.IPM)
                
)



## Set the study years for each population

StudyYears <- list(
  rhpR = 2019:2072,
  ChpR = 2019:2072,
  rhpr = 2019:2072,
  chou = 2019:2072
)



## Extract Tmax for all models\

ModelTmax <- c(length(StudyYears$rhpR),
               length(StudyYears$ChpR),
               length(StudyYears$rhpr),
               length(StudyYears$chou)
)


## Make an empty list for storing estimates for each model

sum.data <- list(
  rhpR = NA,
  ChpR = NA,
  rhpr = NA,
  chou = NA
  
)


## Define a short-cut function for posterior summaries (median, 95% CI)

sam.summary <- function(x){
  unname(quantile(x, probs = c(0.025, 0.5, 0.975)))
}



## For each model, extract relevant measures

ModelIDs <- c('rho x BreedpRob', 'Clutch x BreedpRob', 'rho.inProcess', 'Clutch.onlyPrior')

sum.list <- list()


for(i in 1:length(ModelIDs)){
  
  
  # Set sample matrix
  
  sam.mat <- out.mat[[i]]
  
  
  # Prepare matrices to store results
  NTOT <- matrix(NA, nrow = ModelTmax[i], ncol = 3, dimnames = list(NULL, c('lCI', 'Median', 'uCI')))
  
  
  # Make posterior summary for NTOT in each year
  for(t in 1:ModelTmax[i]){
    
    NTOT[t,] <- sam.summary(sam.mat[,paste0('NTOT[', t, ']')])
    
  }
  
  
  # Organise summaries in a data frame
  model.sum.data <- data.frame(
    ModelID = ModelIDs[i],
    Year = StudyYears[[i]],
    Parameter = rep(c('NTOT'), each = ModelTmax[i])
  )
  
  model.sum.data <- cbind(
    model.sum.data,
    rbind(NTOT))
  
  # Store summary data
  sum.list[[i]] <- model.sum.data
  
} 

sum.list


## Make a combined data frame with all populations
allModel.data <- dplyr::bind_rows(sum.list, .id = "column_label")



## Plot Population sizes

ggplot(allModel.data, aes(x = Year, y = Median)) +
  geom_line(aes(color = ModelID, linetype = Parameter)) +
  geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = ModelID),alpha = 0.2) +
  ggtitle('Total Population Size') + 
  ylab('Estimate') + 
  # facet_wrap(~ModelID, ncol = 1, scales = 'free_y') +
  theme_bw() + theme(legend.position = "right", panel.grid.major.y = element_blank()
                   , panel.grid.minor.y = element_blank(),
                   axis.text.x = element_text(angle = 45, vjust = 0.5),
                   plot.title = element_text(face = 'bold'))




#-----------------------------------------------------------------------------






