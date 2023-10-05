
library(tidyverse)
library(dplyr)

# Count in Non-Breeding Season

# Year Label
NB_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2) + 1 

ChF_NBy <- c(16,26,69,80,60,20,25,24,48,59,11,12)    #Female Chicks Count
# ChF_NBm <- c((16 + 26 + 69 + 80 + 60 + 20), (25 + 24 + 48 + 59 + 11 + 12))

AF_NBy <- c(71,69,134,108,86,33,101,68,71,90,34,13)  #All Female 
# AF_NBm <- c((71 + 69 + 134 + 108 + 86 + 33), (101 + 68 + 71 + 90 + 34 + 13))


# Count in Breeding Season

# Year Label

BN_yr <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3)

JuF_BNy <- c(9,18,33,33,10,13,6,9,13,27,14,25,11,0)   #Female Juvenile count
# JuF_BNm <- c((9 + 18 + 33 + 33 + 10 + 13), (6 + 9 + 13 + 27 + 14 + 25), (11 + 0))

AF_BNy <- c(18,49,79,71,47,44,30,25,38,54,48,97,17,1) #All Female
# AF_BNm <- c((18 + 49 + 79 + 71 + 47 + 44), (30 + 25 + 38 + 54 + 48 + 97), (17 + 1))

#  Male count Data

ChM_NBy <- c(20,39,78,72,52,21,25,26,33,54,13,4)      #Male Chicks Count
# ChM_NBm <- c((20 + 39 + 78 + 72 + 52 + 21), (25 + 26 + 33 + 54 + 13 + 4)) 

JuM_BNy <- c(13,17,28,29,8,8,20,9,10,36,11,19,11,0)   #Male juvenile count
# JuM_BNm <- c((13 + 17 + 28 + 29 + 8 + 8), (20 + 9 + 10 + 36 + 11 + 19), (11 + 0))

M1y_BNy <- c(4,16,11,19,13,18,8,2,9,14,16,28,1,0)     #Male 1 year count in Breeding
# M1y_BNm <- c((4 + 16 + 11 + 19 + 13 + 18), (8 + 2 + 9 + 14 + 16 + 28), (1 + 0))

M1y_NBy <- c(9,6,9,13,5,5,12,19,14,10,5,10)           #Male 1 year count in Non-Breeding
# M1y_NBm <- c((9 + 6 + 9 + 13 + 5 + 5), (12 + 19 + 14 + 10 + 5 + 10)) 

M2y_BNy <- c(1,6,4,6,2,6,4,6,1,8,2,14,0,0)            #Male 2 years count in Breeding
# M2y_BNm <- c((1 + 6 + 4 + 6 + 2 + 6), (4 + 6 + 1 + 8 + 2 + 14), (0 + 0)) 

M2y_NBy <- c(5,1,8,4,0,2,2,1,0,0,0,1)                 #Male 2 years count in Non-Breeding
# M2y_NBm <- c((5 + 1 + 8 + 4 + 0 + 2), (2 + 1 + 0 + 0 + 0 + 1))

M3y_BNy <- c(27,118,164,119,49,57,54,60,76,91,111,83,26,1) #Male 3 years count in Breeding
# M3y_BNm <- c((27 + 118 + 164 + 119 + 49 + 57), (54 + 60 + 76 + 91 + 111 + 83), (26 + 1))

M3y_NBy <- c(75,57,143,139,85,34,71,63,83,82,56,37)   #Male 3 years count in Non-Breeding
# M3y_NBm <- c((75 + 57 + 143 + 139 + 85 + 34), (71 + 63 + 83 + 82 + 56 + 37))



# Start from November 2019 - December 2021, Breeding -> Non-breeding

month <- c('Nov19', 'Dec19', 'Jan20', 'Feb20', 'Mar20' ,'Apr20',
           'May20', 'Jun20',' Jul20', 'Aug20', 'Sep20', 'Oct20',
           'Nov20', 'Dec20', 'Jan21', 'Feb21', 'Mar21' ,'Apr21',
           'May21', 'Jun21',' Jul21', 'Aug21', 'Sep21', 'Oct21', 'Nov21', 'Dec21')

# Female_1 <- c(9,18,33,33,10,13,
#               16,26,69,80,60,20,
#               6,9,13,27,14,25,
#               25,24,48,59,11,12,11,0)
# 
# Female_AD <- c(18,49,79,71,47,44,
#                71,69,134,108,86,33,
#                30,25,38,54,48,97,
#                101,68,71,90,34,13,17,1)
# 
# Male_1 <- c(13,17,28,29,8,8,
#             20,39,78,72,52,21,
#             20,9,10,36,11,19,
#             25,26,33,54,13,4,11,0)
# 
# Male_2 <- c(4,16,11,19,13,18,
#             9,6,9,13,5,5,
#             8,2,9,14,16,28,
#             12,19,14,10,5,10,1,0) 
# 
# 
# Male_3 <- c(1,6,4,6,2,6,
#             5,1,8,4,0,2,
#             4,6,1,8,2,14,
#             2,1,0,0,0,1,0,0)                  
# 
# Male_4 <- c(27,118,164,119,49,57,
#             75,57,143,139,85,34,
#             54,60,76,91,111,83,
#             71,63,83,82,56,37,26,1)                  


Allin <- data.frame(
  PopulationSiz=c(9,18,33,33,10,13,
                  16,26,69,80,60,20,
                  6,9,13,27,14,25,
                  25,24,48,59,11,12,11,0, # Female Chick + Jevenile
                  18,49,79,71,47,44,
                  71,69,134,108,86,33,
                  30,25,38,54,48,97,
                  101,68,71,90,34,13,17,1, # Adult Female
                  13,17,28,29,8,8,
                  20,39,78,72,52,21,
                  20,9,10,36,11,19,
                  25,26,33,54,13,4,11,0, # Male Chick + Jevenile
                  4,16,11,19,13,18,
                  9,6,9,13,5,5,
                  8,2,9,14,16,28,
                  12,19,14,10,5,10,1,0, # Male Yearling
                  1,6,4,6,2,6,
                  5,1,8,4,0,2,
                  4,6,1,8,2,14,
                  2,1,0,0,0,1,0,0,       # Male 2 years
                  27,118,164,119,49,57,
                  75,57,143,139,85,34,
                  54,60,76,91,111,83,
                  71,63,83,82,56,37,26,1 # Male 3 years and older
                  ))                   


viz <- data.frame(
  AgeClass = rep(c('Female Chick', 'Adult Female', 'Male Chick', 'Male Yearling',
                   'Male 2 years', 'Male 3 years'), each = 26),
  Year = c('2019','2019','2020','2020','2020','2020','2020','2020',
           '2020','2020','2020','2020','2020','2020','2021','2021',
           '2021','2021','2021','2021','2021','2021','2021','2021',
           '2021','2021'),
  Month = c('Nov19', 'Dec19', 'Jan20', 'Feb20', 'Mar20' ,'Apr20',
            'May20', 'Jun20', 'Jul20', 'Aug20', 'Sep20', 'Oct20',
            'Nov20', 'Dec20', 'Jan21', 'Feb21', 'Mar21' ,'Apr21',
            'May21', 'Jun21', 'Jul21', 'Aug21', 'Sep21', 'Oct21', 
            'Nov21', 'Dec21'), 
  Season = c('Breeding', 'Breeding', 'Breeding', 'Breeding',
             'Breeding', 'Breeding', 'Non-breeding', 'Non-breeding',
             'Non-breeding', 'Non-breeding', 'Non-breeding', 'Non-breeding',
             'Breeding', 'Breeding', 'Breeding', 'Breeding', 'Breeding', 'Breeding',
             'Non-breeding', 'Non-breeding','Non-breeding', 'Non-breeding',
             'Non-breeding','Non-breeding', 'Breeding', 'Breeding'))                  


pop.survey <- cbind(viz, Allin)

pop.survey$AgeClass <- factor(pop.survey$AgeClass, levels = c('Female Chick', 'Adult Female',
                       'Male Chick', 'Male Yearling','Male 2 years', 'Male 3 years'))

# pop.survey$Year <- factor(pop.survey$Year, levels = c('2019','2019','2020','2020',
#                    '2020','2020','2020','2020','2020','2020','2020','2020','2020','2020',
#                    '2021','2021','2021','2021','2021','2021','2021','2021','2021','2021',
#                    '2021','2021'))

pop.survey$Month <- factor(pop.survey$Month, levels = c('Nov19','Dec19', 
                    'Jan20','Feb20','Mar20','Apr20','May20','Jun20',
                    'Jul20','Aug20','Sep20','Oct20','Nov20','Dec20', 
                    'Jan21','Feb21','Mar21','Apr21','May21','Jun21',
                    'Jul21','Aug21','Sep21', 'Oct21','Nov21','Dec21'))

# pop.survey$Season <- factor(pop.survey$Season, levels = c('Breeding','Breeding',
#                      'Breeding','Breeding','Breeding','Breeding','Non-breeding',
#                      'Non-breeding','Non-breeding','Non-breeding','Non-breeding',
#                      'Non-breeding','Breeding','Breeding','Breeding','Breeding',
#                      'Breeding','Breeding','Non-breeding','Non-breeding',
#                      'Non-breeding','Non-breeding','Non-breeding','Non-breeding',
#                      'Breeding','Breeding'))

ggplot(pop.survey, aes(x = Month, y = PopulationSiz))+
       geom_point(aes(color = AgeClass, group = AgeClass), size = 2)+
       geom_line(aes(color = AgeClass, group = AgeClass), size = 0.8)+
       geom_vline(xintercept = 'May20', linetype = "solid",
            color = "black", size = 0.5) +     
       geom_vline(xintercept = 'Nov20', linetype = "solid",
             color = "black", size = 0.5) +
       geom_vline(xintercept = 'May21', linetype = "solid",
             color = "black", size = 0.5) +    
       geom_vline(xintercept = 'Nov21', linetype = "solid",
             color = "black", size = 0.5) +
       geom_text(aes(x='Jan20', label = "Breeding", y = 175),
             colour="black", size = 3, alpha = 0.5) +
       geom_text(aes(x='Aug20', label = "Non-breeding", y = 175),
             colour="black", size = 3, alpha = 0.5) + 
       geom_text(aes(x='Feb21', label = "Breeding", y = 175),
             colour="black", size = 3, alpha = 0.5) +
       geom_text(aes(x='Aug21', label = "Non-breeding", y = 175),
             colour="black", size = 3, alpha = 0.5) + 
       # geom_rect(aes(xmin = 'Nov19', xmax = 'May20', fill = 'Blue'),
       #       ymin = -Inf, ymax = Inf, alpha = 0.5) +
       # geom_rect(aes(xmin = 'May20', xmax = 'Nov20', fill = Season),
       #       ymin = -Inf, ymax = Inf, alpha = 0.5) +
       # geom_rect(aes(xmin = 'Nov20', xmax = 'May21', fill = Season),
       #       ymin = -Inf, ymax = Inf, alpha = 0.5) +
       # geom_rect(aes(xmin = 'May21', xmax = 'Nov21', fill = Season),
       #       ymin = -Inf, ymax = Inf, alpha = 0.5) +
       # scale_fill_manual(values = c("blue", "red"))
       scale_y_continuous()+
       facet_wrap(~ AgeClass, ncol = 2)+
       ggtitle('Population count in survey period') + 
       ylab("Population Count") + 
       theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                        plot.title = element_text(face = 'bold'))











