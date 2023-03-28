calculate.elasticity <- function(A_orig, dy){
  
  elas <- matrix(0, nrow = nrow(A_orig), ncol = ncol(A_orig), dimnames = dimnames(A_orig))
  
  for(m in 1:nrow(A_orig)){
    for(n in 1:ncol(A_orig)){
      
      if(A_orig[m, n] != 0){ # Only continue if matrix element [m, n] != 0
        
        # Define perturbation matrix  
        A_pert <- A_orig
        
        # Perturb target element in matrix
        A_pert[m, n] <- A_pert[m, n] * (1 + dy) 
        
        # Calculate population growth rate for both matrices
        lam_orig <- as.numeric(eigen(A_orig)$values[1])
        lam_pert <- as.numeric(eigen(A_pert)$values[1])
        
        # Calculate sensitivity of population growth rate to target element
        elas[m, n] <- (lam_pert - lam_orig) / (lam_orig * dy)
      }
    }
  }
  return(elas)
}  
