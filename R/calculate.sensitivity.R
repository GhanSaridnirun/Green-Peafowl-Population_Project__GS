calculate.sensitivity <- function(A_orig, dy){
  
  sens <- matrix(NA, nrow = nrow(A_orig), ncol = ncol(A_orig), dimnames = dimnames(A_orig))
  
  for(j in 1:nrow(A_orig)){
    for(k in 1:ncol(A_orig)){
      
      if(A_orig[j, k] != 0){ # Only continue if matrix element [j, k] != 0
        
        # Define perturbation matrix  
        A_pert <- A_orig
        
        # Perturb target element in matrix
        A_pert[j, k] <- A_pert[j, k] + dy
        
        # Calculate population growth rate for both matrices
        lam_orig <- as.numeric(eigen(A_orig)$values[1])
        lam_pert <- as.numeric(eigen(A_pert)$values[1])
        
        # Calculate sensitivity of population growth rate to target element
        sens[j, k] <- (lam_pert - lam_orig) / dy 
      }
    }
  }
  
  return(sens)
} 