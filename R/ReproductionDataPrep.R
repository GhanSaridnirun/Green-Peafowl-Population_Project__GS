
RepData <- read.csv("R/GP_BreedData.csv")
str(RepData)


# Brood size (Chick)

  RepFCh <- subset(RepData, Symbol == 'F_Chick')
  
    xmax <- length(RepFCh$Brood.size)               # Data length
    
    BroodSize <- as.vector(RepFCh$Brood.size)       # Brood size
    
    Year_BS <- as.vector(RepFCh$index_y)            # Year index Breeding   

      
# Breeding probability
    
  NonBreedRep <- subset(RepData, Season == 'Non-Breeding')
    
  NonBreedRep$Rep <- ifelse(NonBreedRep$Age.Class == 'Single_Female', 0, 1)
  
  table(NonBreedRep$Age.Class, NonBreedRep$Rep)

  Rep <- as.vector(NonBreedRep$Rep)
  
  ymax <- length(Rep)
  
  
# Brood data in breeding season 
  
  BreedRep <- subset(RepData, Season == 'Breeding')
  
  BreedRep$Rep <- ifelse(BreedRep$Age.Class == 'Single_Female', 0, 1)
  
  table(BreedRep$Age.Class, BreedRep$Rep)
  
  BrRep <- as.vector(BreedRep$Rep)
  
  vmax <- length(BrRep)
  

# Clutch Sizes
  
  ClutchSize <- c(3,7,4,4,4,11,3,3,5,5,10)  


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  