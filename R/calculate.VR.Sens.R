calculate.VR.Sens <- function(dy, VR.names, VR.orig){
  
  VR.sens.data <- data.frame()
  for(i in 1:length(VR.names)){ # Loop over each perturbation scenario
    
    # Build perturbed matrix
    A_pert <- make.GPprojMatrix(pRep = VR.pert["pRep", i], 
                                mean.CS = VR.pert["mean.CS", i], 
                                S_C = VR.pert["S_C", i], 
                                gamma = VR.pert["gamma", i], 
                                sF_NB = VR.pert[c("sF_NB_Ch", paste0("sF_NB_", 1:3, "yr")), i], 
                                sF_BN = VR.pert[c("sF_BN_Ju", paste0("sF_BN_", 1:3, "yr")), i], 
                                sM_NB = VR.pert[c("sM_NB_Ch", paste0("sM_NB_", 1:3, "yr")), i], 
                                sM_BN = VR.pert[c("sM_BN_Ju", paste0("sM_BN_", 1:3, "yr")), i],
                                seasonal = FALSE)$A
    
    # Calculate population growth rate for both matrices
    lam_orig <- as.numeric(eigen(A_orig)$values[1])
    lam_pert <- as.numeric(eigen(A_pert)$values[1])
    
    # Calculate sensitivity of population growth rate to target element
    VR.sens <- (lam_pert - lam_orig) / dy
    
    # Assemble results in a dataframe
    data.temp <- data.frame(VitalRate = VR.names[i],
                            Sensitivity = VR.sens)
    
    VR.sens.data <- rbind(VR.sens.data, data.temp)
  }
  return(VR.sens.data)
}