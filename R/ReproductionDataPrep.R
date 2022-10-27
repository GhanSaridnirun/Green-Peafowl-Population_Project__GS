
RepData <- read.csv("R/GP_BreedData.csv")
str(RepData)


# Brood size (Chick)

  RepFCh <- subset(RepData, Symbol == 'F_Chick')
  
    xmax <- length(RepFCh$Brood.size)               # Data length
    
    BroodSize <- as.vector(RepFCh$Brood.size)       # Brood size
    
    Year_BS <- as.vector(RepFCh$index_y)            # Year index Breeding   
  
  

