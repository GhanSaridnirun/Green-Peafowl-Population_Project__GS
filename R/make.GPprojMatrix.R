
make.GPprojMatrix <- function(pRep, mean.CS, S_C, gamma, sF_NB, sF_BN, sM_NB, sM_BN,
                              seasonal = FA.conceptLSE){
  
  if(seasonal){
    
    ## Set up seasonal matrix
    A.concept <- matrix(0, nrow = 16, ncol = 16,
                dimnames = list(c("ChF_NB", "ChM_NB", "JuF_BN", "JuM_BN",
                                  "1yF_NB", "1yM_NB", "1yF_BN", "1yM_BN",
                                  "2yF_NB", "2yM_NB", "2yF_BN", "2yM_BN",
                                  "3yF_NB", "3yM_NB", "3yF_BN", "3yM_BN"),
                        
                                c("ChF_NB", "ChM_NB", "JuF_BN", "JuM_BN",
                                  "1yF_NB", "1yM_NB", "1yF_BN", "1yM_BN",
                                  "2yF_NB", "2yM_NB", "2yF_BN", "2yM_BN",
                                  "3yF_NB", "3yM_NB", "3yF_BN", "3yM_BN")))
    
    A.concept.concept <- A.concept
    
    ## Populate seasonal matrix (numerical)
    # Row index = gets produced
    # Column index = produces
    
    # Female chicks (NB)
    A.concept["ChF_NB", "2yF_BN"] <- pRep * mean.CS * S_C * gamma
    A.concept["ChF_NB", "3yF_BN"] <- pRep * mean.CS * S_C * gamma
    
    # Male chicks (NB)
    A.concept["ChM_NB", "2yF_BN"] <- pRep * mean.CS * S_C * (1-gamma)
    A.concept["ChM_NB", "3yF_BN"] <- pRep * mean.CS * S_C * (1-gamma)
    
    # Female juvenile (BN)
    A.concept["JuF_BN", "ChF_NB"] <- sF_NB[1]
    
    # Male juvenile (BN)
    A.concept["JuM_BN", "ChM_NB"] <- sM_NB[1]
    
    # Female yearling (NB)
    A.concept["1yF_NB", "JuF_BN"] <- sF_BN[1]
    
    # Male yearling (NB)
    A.concept["1yM_NB", "JuM_BN"] <- sM_BN[1]
    
    # Female yearling (BN)
    A.concept["1yF_BN", "1yF_NB"] <- sF_NB[2]
    
    # Male yearling (BN)
    A.concept["1yM_BN", "1yM_NB"] <- sM_NB[2]
    
    # Female 2-yr (NB)
    A.concept["2yF_NB", "1yF_BN"] <- sF_BN[2]
    
    # Male 2-yr (NB)
    A.concept["2yM_NB", "1yM_BN"] <- sM_BN[2]
    
    # Female 2-yr (BN)
    A.concept["2yF_BN", "2yF_NB"] <- sF_NB[3]
    
    # Male 2-yr (BN)
    A.concept["2yM_BN", "2yM_NB"] <- sM_NB[3]
    
    # Female 3-yr (NB)
    A.concept["3yF_NB", "2yF_BN"] <- sF_BN[3]
    A.concept["3yF_NB", "3yF_BN"] <- sF_BN[4]
    
    # Male 3-yr (NB)
    A.concept["3yM_NB", "2yM_BN"] <- sM_BN[3]
    A.concept["3yM_NB", "3yM_BN"] <- sM_BN[4]
    
    # Females 3-yr (BN)
    A.concept["3yF_BN", "3yF_NB"] <- sF_NB[4]
    
    # Males 3-yr (BN)
    A.concept["3yM_BN", "3yM_NB"] <- sM_NB[4]
    
    
    ## Populate seasonal matrix (conceptual)
    
    # Female chicks (NB)
    A.concept.concept["ChF_NB", "2yF_BN"] <- "pRep * mean.CS * S_C * gamma"
    A.concept.concept["ChF_NB", "3yF_BN"] <- "pRep * mean.CS * S_C * gamma"
    
    # Male chicks (NB)
    A.concept.concept["ChM_NB", "2yF_BN"] <- "pRep * mean.CS * S_C * (1-gamma)"
    A.concept.concept["ChM_NB", "3yF_BN"] <- "pRep * mean.CS * S_C * (1-gamma)"
    
    # Female juvenile (BN)
    A.concept.concept["JuF_BN", "ChF_NB"] <- "sF_NB[1]"
    
    # Male juvenile (BN)
    A.concept.concept["JuM_BN", "ChM_NB"] <- "sM_NB[1]"
    
    # Female yearling (NB)
    A.concept.concept["1yF_NB", "JuF_BN"] <- "sF_BN[1]"
    
    # Male yearling (NB)
    A.concept.concept["1yM_NB", "JuM_BN"] <- "sM_BN[1]"
    
    # Female yearling (BN)
    A.concept.concept["1yF_BN", "1yF_NB"] <- "sF_NB[2]"
    
    # Male yearling (BN)
    A.concept.concept["1yM_BN", "1yM_NB"] <- "sM_NB[2]"
    
    # Female 2-yr (NB)
    A.concept.concept["2yF_NB", "1yF_BN"] <- "sF_BN[2]"
    
    # Male 2-yr (NB)
    A.concept.concept["2yM_NB", "1yM_BN"] <- "sM_BN[2]"
    
    # Female 2-yr (BN)
    A.concept.concept["2yF_BN", "2yF_NB"] <- "sF_NB[3]"
    
    # Male 2-yr (BN)
    A.concept.concept["2yM_BN", "2yM_NB"] <- "sM_NB[3]"
    
    # Female 3-yr (NB)
    A.concept.concept["3yF_NB", "2yF_BN"] <- "sF_BN[3]"
    A.concept.concept["3yF_NB", "3yF_BN"] <- "sF_BN[4]"
    
    # Male 3-yr (NB)
    A.concept.concept["3yM_NB", "2yM_BN"] <- "sM_BN[3]"
    A.concept.concept["3yM_NB", "3yM_BN"] <- "sM_BN[4]"
    
    # Females 3-yr (BN)
    A.concept.concept["3yF_BN", "3yF_NB"] <- "sF_NB[4]"
    
    # Males 3-yr (BN)
    A.concept.concept["3yM_BN", "3yM_NB"] <- "sM_NB[4]"
    
    
    
  }else{
    
    ## Set up annual matrix
    A.concept <- matrix(0, nrow = 8, ncol = 8,
                dimnames = list(c("ChF", "ChM",
                                  "1yF", "1yM",
                                  "2yF", "2yM",
                                  "3yF", "3yM"),
                        
                                c("ChF", "ChM",
                                  "1yF", "1yM",
                                  "2yF", "2yM",
                                  "3yF", "3yM")))
    
    A.concept.concept <- A.concept
    
    ## Populate annual matrix (numerical)
    # Row index = gets produced
    # Column index = produces
    
    # Female chicks (NB)
    A.concept["ChF", "1yF"] <- sF_NB[2] * pRep * mean.CS * S_C * gamma
    A.concept["ChF", "2yF"] <- sF_NB[3] * pRep * mean.CS * S_C * gamma
    A.concept["ChF", "3yF"] <- sF_NB[4] * pRep * mean.CS * S_C * gamma
    
    # Male chicks (NB)
    A.concept["ChM", "1yF"] <- sF_NB[2] * pRep * mean.CS * S_C * (1-gamma)
    A.concept["ChM", "2yF"] <- sF_NB[3] * pRep * mean.CS * S_C * (1-gamma)
    A.concept["ChM", "3yF"] <- sF_NB[4] * pRep * mean.CS * S_C * (1-gamma)
    
    # Female yearling (NB)
    A.concept["1yF", "ChF"] <- sF_NB[1] * sF_BN[1]
    
    # Male yearling (NB)
    A.concept["1yM", "ChM"] <- sM_NB[1] * sM_BN[1]
    
    # Female 2-yr (NB)
    A.concept["2yF", "1yF"] <- sF_NB[2] * sF_BN[2]
    
    # Male 2-yr (NB)
    A.concept["2yM", "1yM"] <- sM_NB[2] * sM_BN[2]
    
    # Female 3-yr (NB)
    A.concept["3yF", "2yF"] <- sF_NB[3] * sF_BN[3]
    A.concept["3yF", "3yF"] <- sF_NB[4] * sF_BN[4]
    
    # Male 3-yr (NB)
    A.concept["3yM", "2yM"] <- sM_NB[3] * sM_BN[3]
    A.concept["3yM", "3yM"] <- sM_NB[4] * sM_BN[4]
    
    
    ## Populate annual matrix (conceptual)
    
    # Female chicks (NB)
    A.concept["ChF", "1yF"] <- "sF_NB[2] * pRep * mean.CS * S_C * gamma"
    A.concept["ChF", "2yF"] <- "sF_NB[3] * pRep * mean.CS * S_C * gamma"
    A.concept["ChF", "3yF"] <- "sF_NB[4] * pRep * mean.CS * S_C * gamma"
    
    # Male chicks (NB)
    A.concept["ChM", "1yF"] <- "sF_NB[2] * pRep * mean.CS * S_C * (1-gamma)"
    A.concept["ChM", "2yF"] <- "sF_NB[3] * pRep * mean.CS * S_C * (1-gamma)"
    A.concept["ChM", "3yF"] <- "sF_NB[4] * pRep * mean.CS * S_C * (1-gamma)"
    
    # Female yearling (NB)
    A.concept["1yF", "ChF"] <- "sF_NB[1] * sF_BN[1]"
    
    # Male yearling (NB)
    A.concept["1yM", "ChM"] <-" sM_NB[1] * sM_BN[1]"
    
    # Female 2-yr (NB)
    A.concept["2yF", "1yF"] <- "sF_NB[2] * sF_BN[2]"
    
    # Male 2-yr (NB)
    A.concept["2yM", "1yM"] <- "sM_NB[2] * sM_BN[2]"
    
    # Female 3-yr (NB)
    A.concept["3yF", "2yF"] <- "sF_NB[3] * sF_BN[3]"
    A.concept["3yF", "3yF"] <-" sF_NB[4] * sF_BN[4]"
    
    # Male 3-yr (NB)
    A.concept["3yM", "2yM"] <- "sM_NB[3] * sM_BN[3]"
    A.concept["3yM", "3yM"] <-" sM_NB[4] * sM_BN[4]"
    
  }
  
  
  return(list(A = A, A.concept = A.concept))
  
}


