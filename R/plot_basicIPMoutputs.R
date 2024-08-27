plot_basicIPMoutputs <- function(mcmc.out, GP.IPMconstants, GP.IPMdata, estimate.rho){
  
  Tmax <- GP.IPMconstants$ny.data
  
  ## Convert samples to matrix
  out.mat <- as.matrix(mcmc.out)
  
  ## Calculate total population sizes
  for(t in 1:Tmax){
    Ntot_Breed <- rowSums(out.mat[, c(paste0("NBreedF[", 1:4, ", ", t, "]"), paste0("NBreedF[", 1:4, ", ", t, "]"))])
    Ntot_NonBreed <- rowSums(out.mat[, c(paste0("NNonF[", 1:4, ", ", t, "]"), paste0("NNonF[", 1:4, ", ", t, "]"))])
    
    out.mat <- cbind(out.mat, Ntot_Breed, Ntot_NonBreed)
    colnames(out.mat)[ncol(out.mat)-c(1, 0)] <- c(paste0("Ntot_Breed[", t, "]"), paste0("Ntot_NonBreed[", t, "]"))
  }
  
  ## Calculate population proportions
  for(t in 1:Tmax){
    
    Nprops_Breed <- out.mat[, c(paste0("NBreedF[", 1:4, ", ", t, "]"), paste0("NBreedM[", 1:4, ", ", t, "]"))]/out.mat[, paste0("Ntot_Breed[", t, "]")]
    colnames(Nprops_Breed) <- paste0("prop_", colnames(Nprops_Breed))
    
    Nprops_NonBreed <- out.mat[, c(paste0("NNonF[", 1:4, ", ", t, "]"), paste0("NNonM[", 1:4, ", ", t, "]"))]/out.mat[, paste0("Ntot_NonBreed[", t, "]")]
    colnames(Nprops_NonBreed) <- paste0("prop_", colnames(Nprops_NonBreed))
    
    out.mat <- cbind(out.mat, Nprops_Breed, Nprops_NonBreed)
  }
  
  ## List relevant parameters
  
  # Survival probabilities
  survParams <- c(paste0("sF_NB[", 1:4, "]"), paste0("sF_BN[", 1, "]"))
  
  # Reproduction parameters
  repParams <- c("pRep", "mean.CS", "mean.rho", "mean.S_C")
  
  # Total population size
  NtotParams <- c(paste0("Ntot_Breed[", 1:Tmax, "]"), paste0("Ntot_NonBreed[", 2:Tmax, "]"))
  
  # Population proportions
  NpropParams <- colnames(out.mat)[which(stringr::str_detect(colnames(out.mat), "prop"))]
  
  
  ## Arrange results in data frames
  
  # Vital rates
  VR_results <- reshape2::melt(out.mat[, c(survParams, repParams)]) %>%
    dplyr::mutate(Parameter = dplyr::case_when(stringr::str_detect(Var2, "sF") ~ "Seasonal survival",
                                               Var2 == "pRep" ~ "Breeding probability",
                                               Var2 == "mean.CS" ~ "Clutch size",
                                               Var2 == "mean.rho" ~ "Brood size",
                                               Var2 == "mean.S_C" ~ "Survival to fledging"),
                  AgeClassIdx = stringr::str_extract_all(Var2, pattern = "\\d+", simplify = TRUE),
                  Season = dplyr::case_when(stringr::str_detect(Var2, "NB") ~ "NonBreed",
                                            stringr::str_detect(Var2, "BN") ~ "Breed")) %>%
    dplyr::mutate(AgeClass = dplyr::case_when(AgeClassIdx == 1 & Season == "NonBreed" ~ "Chick",
                                              AgeClassIdx == 1 & Season == "Breed" ~ "Juvenile",
                                              AgeClassIdx %in% c(2, 3) ~ "1-2 years",
                                              AgeClassIdx == 4 ~ "3+ years")) %>%
    dplyr::select(Parameter, Season, AgeClass, value)
  
  VR_results$AgeClass <- factor(VR_results$AgeClass, levels = c("Chick", "Juvenile", "1-2 years", "3+ years"))
  
  
  # Total population size
  Ntot_results <- reshape2::melt(out.mat[, NtotParams]) %>%
    dplyr::mutate(Season = dplyr::case_when(stringr::str_detect(Var2, "NonBreed") ~ "NonBreed",
                                            stringr::str_detect(Var2, "Breed") ~ "Breed"),
                  Year = as.numeric(stringr::str_extract_all(Var2, pattern = "\\d+", simplify = TRUE)) + 2019 - 1) %>%
    dplyr::select(Year, Season, value)
  
  
  # Population proportions
  Idx.data <- data.frame(stringr::str_extract_all(NpropParams, pattern = "\\d+", simplify = TRUE)) %>%
    dplyr::rename(AgeClassIdx = X1, YearIdx = X2) %>%
    dplyr::mutate(Var2 = NpropParams)
  
  Nprop_results <- reshape2::melt(out.mat[, NpropParams]) %>%
    dplyr::left_join(Idx.data, by = "Var2") %>%
    plyr::mutate(Season = dplyr::case_when(stringr::str_detect(Var2, "Non") ~ "NonBreed",
                                           stringr::str_detect(Var2, "Breed") ~ "Breed"),
                 Sex = dplyr::case_when(stringr::str_detect(Var2, "edF") | stringr::str_detect(Var2, "onF") ~ "Female",
                                        stringr::str_detect(Var2, "edM") | stringr::str_detect(Var2, "onM") ~ "Male"),
                 AgeClass = dplyr::case_when(AgeClassIdx == 1 ~ "Chick/Juvenile",
                                             AgeClassIdx == 2 ~ "Yearling",
                                             AgeClassIdx == 3 ~ "2 years",
                                             AgeClassIdx == 4 ~ "3+ years"),
                 Year = as.numeric(YearIdx) + 2019 - 1) %>%
    dplyr::select(Year, Season, AgeClass, Sex, value)
  
  Nprop_results$AgeClass <- factor(Nprop_results$AgeClass, levels = c("Chick/Juvenile", "Yearling", "2 years", "3+ years"))
  
  
  ## Simulate vital rate priors
  nsamples <- nrow(out.mat)
  prior_mean.CS <- runif(nsamples, 3, 11)
  prior_mean.rho <- runif(nsamples, 0, 10)
  prior_mean.S_C <- runif(nsamples, 0, 1)
  prior_rRep <- runif(nsamples, 0, 1)
  
  prior_s_yr_ad <- truncnorm::rtruncnorm(nsamples, a = 0, b = 1, mean = metaSurv.mean, sd = metaSurv.se)
  prior_s_yr_sa <- truncnorm::rtruncnorm(nsamples, a = 0, b = 1, mean = metaSurv.mean - saSurv.diff, sd = metaSurv.se)
  prior_s <- runif(nsamples, 0, 1)
  
  prior_data <- data.frame(
    value = c(prior_s, prior_s, sqrt(prior_s_yr_sa), sqrt(prior_s_yr_ad),
              prior_rRep, prior_mean.CS, prior_mean.S_C, prior_mean.rho),
    Parameter = rep(c(rep("Seasonal survival", 4), "Breeding probability", "Clutch size", "Survival to fledging", "Brood size"), each = nsamples),
    AgeClass = rep(c("Chick", "Juvenile", "1-2 years", "3+ years", rep(NA, 4)), each = nsamples)
  )  
  
  prior_data$AgeClass <- factor(prior_data$AgeClass, levels = c("Chick", "Juvenile", "1-2 years", "3+ years"))
  
  ## Assemble count data
  NB_count <- BN_count <- rep(NA, Tmax)
  for(tt in 1:Tmax){
    
    if(t %in% GP.IPMconstants$NB_yr){
      NB_count[t] <- 
        sum(GP.IPMdata$ChF_NB[which(GP.IPMconstants$NB_yr == t)]) + # Female chicks
        sum(GP.IPMdata$AF_NB[which(GP.IPMconstants$NB_yr == t)]) + # Adult females
        sum(colSums(GP.IPMdata$M_NB)[which(GP.IPMconstants$NB_yr == t)]) # All males
    }
    
    BN_count[t] <- 
      sum(GP.IPMdata$JuF_BN[which(GP.IPMconstants$BN_yr == t)]) + # Female juveniles
      sum(GP.IPMdata$AF_BN[which(GP.IPMconstants$BN_yr == t)]) + # Adult females
      sum(colSums(GP.IPMdata$M_BN)[which(GP.IPMconstants$BN_yr == t)]) # All males
  }
  
  NB_count <- GP.IPMdata$ChF_NB + GP.IPMdata$AF_NB + colSums(GP.IPMdata$M_NB)
  BN_count <- GP.IPMdata$JuF_BN + GP.IPMdata$AF_BN + colSums(GP.IPMdata$M_BN)
  
  count.data <- data.frame(
    Year = c(GP.IPMconstants$BN_yr, GP.IPMconstants$NB_yr) + 2019 - 1,
    Season = c(rep("Breed", length(BN_count)), rep("NonBreed", length(NB_count))), 
    Count = c(BN_count, NB_count)
  )
  
  ## Plot results
  
  # Vital rates
  survPanel <- ggplot(subset(VR_results, Parameter == "Seasonal survival"), aes(x = value)) + 
    geom_density(color = "hotpink", fill = "hotpink", alpha = 0.5) + 
    geom_density(data = subset(prior_data, Parameter == "Seasonal survival"), aes(x = value), color = "#A40167", linetype = "dashed") + 
    facet_wrap(~AgeClass, nrow = 1, scales = "free") + 
    ggtitle("A) Survival probabilities") + 
    xlab("Density") + ylab("Estimate") +
    theme_classic()
  
  if(estimate.rho){
    prior_data <- subset(prior_data, Parameter != "Survival to fledging")
  }else{
    prior_data <- subset(prior_data, Parameter != "Brood size")
  }
  
  repPanel <-  ggplot(subset(VR_results, Parameter != "Seasonal survival"), aes(x = value)) + 
    geom_density(color = "hotpink", fill = "hotpink", alpha = 0.5) + 
    geom_density(data = subset(prior_data, Parameter != "Seasonal survival"), aes(x = value), color = "#A40167", linetype = "dashed") + 
    facet_wrap(~Parameter, nrow = 1, scales = "free") + 
    ggtitle("B) Reproductive parameters") + 
    xlab("Density") + ylab("Estimate") +
    theme_classic()
  
  pdf("IPMResults_VitalRates.pdf", width = 8, height = 5)
    gridExtra::grid.arrange(survPanel, repPanel, ncol = 1)
  dev.off()
  
  
  # Population size
  Ntot_panel <- ggplot(Ntot_results, aes(y = as.factor(Year), x = value)) + 
    ggridges::geom_density_ridges(aes(color = Season, fill = Season), alpha = 0.5) +
    geom_point(data = count.data, aes(y = as.factor(Year), x = Count*0.9, color = Season)) + 
    scale_color_manual(values = c("#A40167", "#E58434")) + 
    scale_fill_manual(values = c("#A40167", "#E58434")) + 
    #scale_y_discrete(limits = unique(rev(as.factor(Ntot_results$Year)))) +
    xlab("Population size") + ylab("Year") +
    ggtitle("A) Population size") + 
    theme_classic()
  
  Nprop_panel <-  ggplot(Nprop_results, aes(y = as.factor(Year), x = value)) + 
    ggridges::geom_density_ridges(aes(fill = AgeClass, linetype = Sex), alpha = 0.5, bandwidth = 0.01, scale = 1) +
    scale_fill_viridis_d(option = "mako") + 
    #scale_y_discrete(limits = unique(rev(as.factor(Ntot_results$Year)))) +
    xlab("Prop. of population") + ylab("Year") + 
    xlim(0, 0.65) + 
    facet_wrap(~Season, ncol = 2) + 
    ggtitle("B) Population structure") + 
    theme_classic()
  
  pdf("IPMResults_PopParams.pdf", width = 8, height = 8)
   gridExtra::grid.arrange(Ntot_panel, Nprop_panel, ncol = 1)
  dev.off()
  
}
