pertSurv.nimble <- nimbleFunction(
  
  run = function(Surv = double(0),
                 pertFac_t = double(1)){
    
    # Calculate perturbed survival
    Surv_t_pert <- Surv * pertFac_t
    
    # Replace values above 1 and below 0
    for(t in 1:length(Surv_t_pert)){
      
      if(Surv_t_pert[t] > 1){
        Surv_t_pert[t] <- 1
      }
      
      if(Surv_t_pert[t] < 0){
        Surv_t_pert[t] <- 0
      }
    }
    
    # Return estimate
    return(Surv_t_pert)
    returnType(double(1))
  
})


## Test uncompiled function
#pertSurv.nimble(Surv_t = 0.2, pertFac_t = c(1, 1, 0.9, -1.2, -2))
#pertSurv.nimble(Surv_t = 0.9, pertFac_t = c(1, 1, 0.9, 1.2, 2))

## Compile function
#pertSurv.nimble.c <- compileNimble(pertSurv.nimble)

## Test compiled function
#pertSurv.nimble.c(Surv_t = 0.2, pertFac_t = c(1, 1, 0.9, -1.2, -2))
#pertSurv.nimble.c(Surv_t = 0.9, pertFac_t = c(1, 1, 0.9, 1.2, 2))
