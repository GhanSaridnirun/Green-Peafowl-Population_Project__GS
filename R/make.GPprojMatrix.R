
make.GPprojMatrix <- function(pRep, mean.CS, S_C, gamma, sF_NB, sF_BN, sM_NB, sM_BN,
                              seasonal = FALSE){
  
  if(seasonal){
    
    ## Set up seasonal matrices
    A_BN <- A_NB <- matrix(0, nrow = 16, ncol = 16,
                           dimnames = list(c("ChF_NB", "ChM_NB", "JuF_BN", "JuM_BN",
                                             "1yF_NB", "1yM_NB", "1yF_BN", "1yM_BN",
                                             "2yF_NB", "2yM_NB", "2yF_BN", "2yM_BN",
                                             "3yF_NB", "3yM_NB", "3yF_BN", "3yM_BN"),
                                           
                                           c("ChF_NB", "ChM_NB", "JuF_BN", "JuM_BN",
                                             "1yF_NB", "1yM_NB", "1yF_BN", "1yM_BN",
                                             "2yF_NB", "2yM_NB", "2yF_BN", "2yM_BN",
                                             "3yF_NB", "3yM_NB", "3yF_BN", "3yM_BN")))
    
    A_BN.concept <- A_BN
    A_NB.concept <- A_NB
    
    ## Populate B->N seasonal matrix (numerical)
    # Row index = gets produced
    # Column index = produces
    
    # Female chicks (NB)
    A_BN["ChF_NB", "2yF_BN"] <- pRep * mean.CS * S_C * gamma
    A_BN["ChF_NB", "3yF_BN"] <- pRep * mean.CS * S_C * gamma
    
    # Male chicks (NB)
    A_BN["ChM_NB", "2yF_BN"] <- pRep * mean.CS * S_C * (1-gamma)
    A_BN["ChM_NB", "3yF_BN"] <- pRep * mean.CS * S_C * (1-gamma)
    
    # Female yearling (NB)
    A_BN["1yF_NB", "JuF_BN"] <- sF_BN[1]
    
    # Male yearling (NB)
    A_BN["1yM_NB", "JuM_BN"] <- sM_BN[1]
    
    # Female 2-yr (NB)
    A_BN["2yF_NB", "1yF_BN"] <- sF_BN[2]
    
    # Male 2-yr (NB)
    A_BN["2yM_NB", "1yM_BN"] <- sM_BN[2]
    
    # Female 3-yr (NB)
    A_BN["3yF_NB", "2yF_BN"] <- sF_BN[3]
    A_BN["3yF_NB", "3yF_BN"] <- sF_BN[4]
    
    # Male 3-yr (NB)
    A_BN["3yM_NB", "2yM_BN"] <- sM_BN[3]
    A_BN["3yM_NB", "3yM_BN"] <- sM_BN[4]
 
    
    ## Populate B->N seasonal matrix (conceptual)
    
    # Female chicks (NB)
    A_BN.concept["ChF_NB", "2yF_BN"] <- "pRep * mean.CS * S_C * gamma"
    A_BN.concept["ChF_NB", "3yF_BN"] <- "pRep * mean.CS * S_C * gamma"
    
    # Male chicks (NB)
    A_BN.concept["ChM_NB", "2yF_BN"] <- "pRep * mean.CS * S_C * (1-gamma)"
    A_BN.concept["ChM_NB", "3yF_BN"] <- "pRep * mean.CS * S_C * (1-gamma)"
    
    # Female yearling (NB)
    A_BN.concept["1yF_NB", "JuF_BN"] <- "sF_BN[1]"
    
    # Male yearling (NB)
    A_BN.concept["1yM_NB", "JuM_BN"] <- "sM_BN[1]"
    
    # Female 2-yr (NB)
    A_BN.concept["2yF_NB", "1yF_BN"] <- "sF_BN[2]"
    
    # Male 2-yr (NB)
    A_BN.concept["2yM_NB", "1yM_BN"] <- "sM_BN[2]"
    
    # Female 3-yr (NB)
    A_BN.concept["3yF_NB", "2yF_BN"] <- "sF_BN[3]"
    A_BN.concept["3yF_NB", "3yF_BN"] <- "sF_BN[4]"
    
    # Male 3-yr (NB)
    A_BN.concept["3yM_NB", "2yM_BN"] <- "sM_BN[3]"
    A_BN.concept["3yM_NB", "3yM_BN"] <- "sM_BN[4]"
    
    
    ## Populate N->B seasonal matrix (numerical)
    # Row index = gets produced
    # Column index = produces
    
    # Female juvenile (BN)
    A_NB["JuF_BN", "ChF_NB"] <- sF_NB[1]
    
    # Male juvenile (BN)
    A_NB["JuM_BN", "ChM_NB"] <- sM_NB[1]
    
    # Female yearling (BN)
    A_NB["1yF_BN", "1yF_NB"] <- sF_NB[2]
    
    # Male yearling (BN)
    A_NB["1yM_BN", "1yM_NB"] <- sM_NB[2]
    
    # Female 2-yr (BN)
    A_NB["2yF_BN", "2yF_NB"] <- sF_NB[3]
    
    # Male 2-yr (BN)
    A_NB["2yM_BN", "2yM_NB"] <- sM_NB[3]
    
    # Females 3-yr (BN)
    A_NB["3yF_BN", "3yF_NB"] <- sF_NB[4]
    
    # Males 3-yr (BN)
    A_NB["3yM_BN", "3yM_NB"] <- sM_NB[4]
    
    
    ## Populate N->B seasonal matrix (conceptual)
    
    # Female juvenile (BN)
    A_NB.concept["JuF_BN", "ChF_NB"] <- sF_NB[1]
    
    # Male juvenile (BN)
    A_NB.concept["JuM_BN", "ChM_NB"] <- sM_NB[1]
    
    # Female yearling (BN)
    A_NB.concept["1yF_BN", "1yF_NB"] <- sF_NB[2]
    
    # Male yearling (BN)
    A_NB.concept["1yM_BN", "1yM_NB"] <- sM_NB[2]
    
    # Female 2-yr (BN)
    A_NB.concept["2yF_BN", "2yF_NB"] <- sF_NB[3]
    
    # Male 2-yr (BN)
    A_NB.concept["2yM_BN", "2yM_NB"] <- sM_NB[3]
    
    # Females 3-yr (BN)
    A_NB.concept["3yF_BN", "3yF_NB"] <- sF_NB[4]
    
    # Males 3-yr (BN)
    A_NB.concept["3yM_BN", "3yM_NB"] <- sM_NB[4]
    
    
    ## Assemble matrix list
    projMatrix.list <- list(A_BtoN = A_BN, A_BtoN.concept = A_BN.concept,
                            A_NtoB = A_NB, A_NtoB.concept = A_NB.concept)
    
  }else{
    
    ## Set up annual matrix
    A <- matrix(0, nrow = 8, ncol = 8,
                dimnames = list(c("ChF", "ChM",
                                  "1yF", "1yM",
                                  "2yF", "2yM",
                                  "3yF", "3yM"),
                        
                                c("ChF", "ChM",
                                  "1yF", "1yM",
                                  "2yF", "2yM",
                                  "3yF", "3yM")))
    
    A.concept <- A
    
    ## Populate annual matrix (numerical)
    # Row index = gets produced
    # Column index = produces
    
    # Female chicks (NB)
    A["ChF", "2yF"] <- sF_NB[3] * pRep * mean.CS * S_C * gamma
    A["ChF", "3yF"] <- sF_NB[4] * pRep * mean.CS * S_C * gamma
    
    # Male chicks (NB)
    A["ChM", "2yF"] <- sF_NB[3] * pRep * mean.CS * S_C * (1-gamma)
    A["ChM", "3yF"] <- sF_NB[4] * pRep * mean.CS * S_C * (1-gamma)
    
    # Female yearling (NB)
    A["1yF", "ChF"] <- sF_NB[1] * sF_BN[1]
    
    # Male yearling (NB)
    A["1yM", "ChM"] <- sM_NB[1] * sM_BN[1]
    
    # Female 2-yr (NB)
    A["2yF", "1yF"] <- sF_NB[2] * sF_BN[2]
    
    # Male 2-yr (NB)
    A["2yM", "1yM"] <- sM_NB[2] * sM_BN[2]
    
    # Female 3-yr (NB)
    A["3yF", "2yF"] <- sF_NB[3] * sF_BN[3]
    A["3yF", "3yF"] <- sF_NB[4] * sF_BN[4]
    
    # Male 3-yr (NB)
    A["3yM", "2yM"] <- sM_NB[3] * sM_BN[3]
    A["3yM", "3yM"] <- sM_NB[4] * sM_BN[4]
    
    
    ## Populate annual matrix (conceptual)
    
    # Female chicks (NB)
    A.concept["ChF", "2yF"] <- "sF_NB[3] * pRep * mean.CS * S_C * gamma"
    A.concept["ChF", "3yF"] <- "sF_NB[4] * pRep * mean.CS * S_C * gamma"
    
    # Male chicks (NB)
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
    
    
    ## Assemble matrix list
    projMatrix.list <- list(A = A, A.concept = A.concept)
  }
  
  
  return(projMatrix.list)
  
}


