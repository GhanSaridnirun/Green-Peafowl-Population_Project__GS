
# Meleagris gallopavo (Wild Turkey) 

## Reported value for Female Survival

MgF01 <- c(0.557-0.172,0.557+0.172) # Average of MgF01 - MgF14            #Ch10
MgF02 <- c(0.30,0.76)                                                     #Ch11
MgF03 <- c(0.842-0.032,0.842+0.032) # Average of MgF15 - MgF21            #Ch07
MgF04 <- c(0.37-0.15,0.37+0.15)                                           #Ch26
MgF05 <- c(0.58-0.11,0.58+0.11)                                           #Ch23
MgF06 <- c(0.649-0.053,0.649+0.053)                                       #Ch21
MgF07 <- c(0.68-0.02,0.68+0.02)                                           #Ch20


all.MgF <- rbind(MgF01,MgF02,MgF03,MgF04,MgF05,MgF06,MgF07)

mean(all.MgF) 
median(all.MgF) 
min(all.MgF) 
max(all.MgF) 

quantile(all.MgF, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))



## Reported value for Male Survival


MgM01 <- c(0.775-0.068,0.775+0.068)  # Average of MgM01 - MgF06
MgM02 <- c(0.49-0.04,0.49+0.04)                                         #Ch20

all.MgM <- rbind(MgM01,MgM02)

mean(all.MgM) 
median(all.MgM) 
min(all.MgM) 
max(all.MgM)

quantile(all.MgM, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

##------------------------------------------------------------


# Tetrao urogallus (Western Capercaillie) 


TuBJu01 <- c(0.32,0.66) # Both Sex Juvenile                        #Ch16
TuBAd01 <- c(0.60,0.81) # Both Sex Adult                           #Ch16
TuF01 <- c(0.62,0.75)   # Female                                   #Ch18

all.Tu <- rbind(TuBJu01,TuBAd01,TuF01)

mean(all.Tu) 
median(all.Tu) 
min(all.Tu) 
max(all.Tu) 

quantile(all.Tu, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))


# (Centrocercus urophasianus) Greater Sage-Grouse


CuFY01 <- c(0.61,0.69) # Female Yearling                          #Ch14
CuFA01 <- c(0.54,0.61) # Female Adult                             #Ch14

CuA01 <- c(0.090,0.398) # Both Sex Yearling+Adult 2013-2014       #Ch13
CuA02 <- c(0.163,0.558) # Both Sex Yearling+Adult 2014-2015       #Ch13

all.Cu <- rbind(CuFY01,CuFA01,CuA01,CuA02)

mean(all.Cu) 
median(all.Cu) 
min(all.Cu) 
max(all.Cu) 

quantile(all.Cu, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))



# Visualization 


par(mfrow = c(1,2), mar=c(5.1,3.0,3.0,1))
plot(MgF01, rep(0.1, 2), type = 'l', xlim = c(0, 1.0), ylim = c(0, 0.8), xlab = 'Reported values', yaxt = 'n', ylab = '', lwd = 1.5)
title('a) Annual Female Wild Turkey survival', adj = 0)
text(max(MgF01)+0.07, 0.1,'MgF01', col ='grey')
lines(MgF02, rep(0.2, 2), lwd = 1.5)
text(max(MgF02)+0.07, 0.2,'MgF02', col ='grey')
lines(MgF03, rep(0.3, 2), lwd = 1.5)
text(max(MgF03)+0.07, 0.3,'MgF03', col ='grey')
lines(MgF04, rep(0.4, 2), lwd = 1.5)
text(max(MgF04)+0.07, 0.4,'MgF04', col ='grey')
lines(MgF05, rep(0.5, 2), lwd = 1.5)
text(max(MgF05)+0.07, 0.5,'MgF05', col ='grey')
lines(MgF06, rep(0.6, 2), lwd = 1.5)
text(max(MgF06)+0.07, 0.6,'MgF06', col ='grey')
lines(MgF07, rep(0.7, 2), lwd = 1.5)
text(max(MgF07)+0.07, 0.7,'MgF07', col ='grey')

abline(v = quantile(all.MgF, probs = c(0.25, 0.75)), lty = 2, col = '#00C0E2')
abline(v = mean(all.MgF), col = '#ca0020', lty = 2)


plot(MgM01, rep(0.1, 2), type = 'l', xlim = c(0, 0.9), ylim = c(0, 0.8), xlab = 'Reported values', yaxt = 'n', ylab = '', lwd = 1.5)
title('b) Annual Male Wild Turkey survival', adj = 0)
text(max(MgM01)+0.07, 0.1,'MgM01', col ='grey')
lines(MgM02, rep(0.2, 2), lwd = 1.5)
text(max(MgM02)+0.07, 0.2,'MgM02', col ='grey')

abline(v = quantile(all.MgM, probs = c(0.25, 0.75)), lty = 2, col = '#00C0E2')
abline(v = mean(all.MgM), col = '#ca0020', lty = 2)


plot(TuBJu01, rep(0.1, 2), type = 'l', xlim = c(0, 0.9), ylim = c(0, 0.8), xlab = 'Reported values', yaxt = 'n', ylab = '', lwd = 1.5)
title('c) Annual Capercaillie survival', adj = 0)
text(max(TuBJu01)+0.07, 0.1,'TuBJu01', col ='grey')
lines(TuBAd01, rep(0.2, 2), lwd = 1.5)
text(max(TuBAd01)+0.07, 0.2,'TuBAd01', col ='grey')
lines(TuF01, rep(0.3, 2), lwd = 1.5)
text(max(TuF01)+0.07, 0.3,'TuF01', col ='grey')

abline(v = quantile(all.Tu, probs = c(0.25, 0.75)), lty = 2, col = '#00C0E2')
abline(v = mean(all.Tu), col = '#ca0020', lty = 2)


plot(CuFY01, rep(0.1, 2), type = 'l', xlim = c(0, 0.9), ylim = c(0, 0.8), xlab = 'Reported values', yaxt = 'n', ylab = '', lwd = 1.5)
title('d) Annual Greater Sage-Grouse survival', adj = 0)
text(max(CuFY01)+0.07, 0.1,'CuFY01', col ='grey')
lines(CuFA01, rep(0.2, 2), lwd = 1.5)
text(max(CuFA01)+0.07, 0.2,'CuFA01', col ='grey')
lines(CuA01, rep(0.3, 2), lwd = 1.5)
text(max(CuA01)+0.07, 0.3,'CuA01', col ='grey')
lines(CuA02, rep(0.4, 2), lwd = 1.5)
text(max(CuA02)+0.07, 0.4,'CuA02', col ='grey')

abline(v = quantile(all.Cu, probs = c(0.25, 0.75)), lty = 2, col = '#00C0E2')
abline(v = mean(all.Cu), col = '#ca0020', lty = 2)





##------------------------------------------------------------



## All Data Bundle Meleagris gallopavo (Wild Turkey) Survival report


## Reported value for post-hatch - Chick Survival
# MgCh01 <- c(0.69)                                                       #Ch04
MgCh02 <- c(1-0.70,1-0.42) # 0-4 weeks                                    #Ch10
MgCh03 <- c(0.52-0.14,0.52+0.14) # Telemetry                              #Ch01
MgCh04 <- c(0.40-0.15,0.40+0.15) # Flush                                  #Ch01
MgCh05 <- c(0.11,0.67) # 0-16 days                                        #Ch08
MgCh06 <- c(1-0.75,1-0.50) # 2 weeks                                      #Ch11
# MgCh07 <- c(0.64) # 8 weeks Adutl rear                                  #Ch11
# MgCh08 <- c(0.57) # 8 weeks Sub-Adutl rear                              #Ch11
MgCh09 <- c(0.11,0.67) # 0-16 days                                        #Ch03
MgCh10 <- c(0.12,0.52) # 0-2 weeks                                        #Ch02
MgCh11 <- c(0.379-0.111,0.379+0.111) # 0-4 weeks  2010                    #Ch06
MgCh12 <- c(0.352-0.063,0.352+0.063) # 0-4 weeks  2010                    #Ch06
MgCh12.5 <- c(0.366-0.087,0.366+0.087) # Average of MgCh11 & MgCh12
MgCh13 <- c(0.316-0.118,0.316+0.118)  # 0-4 weeks  Forest                 #Ch07  
MgCh14 <- c(0.421-0.043,0.421+0.043)  # 0-4 weeks  Open                   #Ch07
MgCh14.5 <- c(0.369-0.081,0.369+0.081) # Average of MgCh13 & MgCh14
MgCh15 <- c(0.33-0.05,0.33+0.05)  # Adult rear                            #Ch09
MgCh16 <- c(0.16-0.05,0.16+0.05)  # Yearling rear                         #Ch09


## Reported value for Juvenile Survival

### Over all from Ch07
MgJu01 <- c(0.788-0.129,0.788+0.129) # 28 days - 9 months Forest          #Ch07
MgJu02 <- c(0.931-0.028,0.931+0.028) # 28 days - 9 months Open            #Ch07
MgJu03 <- c(0.860-0.079,0.860+0.079) # Average of MgJu01 & MgJu02


## Reported value for Female Juvenile Survival

MgJuF01 <- c(0.667-0.064,0.667+0.064) # 28 days to 9 months Forest-Spring(1 Mar–31 May: 3 months)      #Ch07  
MgJuF02 <- c(0.535-0.157,0.535+0.157) # 28 days to 9 months Forest-Summer(1 Jun–31 Aug: 3 months)      #Ch07
MgJuF03 <- c(0.775-0.030,0.775+0.030) # 28 days to 9 months Forest-Fall-Winter(1 Sep–28 Feb: 6 months) #Ch07
MgJuF04 <- c(0.886-0.035,0.886+0.035) # 28 days to 9 months Open-Spring(1 Mar–31 May: 3 months)        #Ch07   
MgJuF05 <- c(0.720-0.112,0.720+0.112) # 28 days to 9 months Open-Summer(1 Jun–31 Aug: 3 months)        #Ch07
MgJuF06 <- c(0.829-0.171,0.829+0.171) # 28 days to 9 months Open-Fall-Winter(1 Sep–28 Feb: 6 months)   #Ch07
MgJuF06.5 <- c(0.735-0.095,0.735+0.095) # Average of MgJu01 - MgJu06
MgJuF07 <- c(0.98-0.01,0.98+0.01)                                         #Ch20


## Reported value for Male Juvenile Survival
MgJuM01 <- c(0.844-0.032,0.844+0.032) # 28 days to 9 months Forest-Spring(1 Mar–31 May: 3 months)      #Ch07   
MgJuM02 <- c(0.885-0.026,0.885+0.026) # 28 days to 9 months Forest-Summer(1 Jun–31 Aug: 3 months)      #Ch07 
MgJuM03 <- c(0.917-0.021,0.917+0.021) # 28 days to 9 months Forest-Fall-Winter(1 Sep–28 Feb: 6 months) #Ch07
MgJuM04 <- c(0.653-0.027,0.653+0.027) # 28 days to 9 months Open-Spring(1 Mar–31 May: 3 months)        #Ch07   
MgJuM05 <- c(0.937-0.019,0.937+0.019) # 28 days to 9 months Open-Summer(1 Jun–31 Aug: 3 months)        #Ch07 
MgJuM06 <- c(0.872-0.021,0.872+0.021) # 28 days to 9 months Open-Fall-Winter(1 Sep–28 Feb: 6 months)   #Ch07
MgJuM06.5 <- c(0.827-0.056,0.827+0.056) # Average of MgJuM01 - MgJuM06
MgJuM07 <- c(0.96-0.14,0.96+0.14)                                         #Ch20


## Reported value for Female Survival

MgF01 <- c(0.528-0.149,0.528+0.149)                                       #Ch10
MgF02 <- c(0.549-0.482,0.549+0.482)                                       #Ch10
MgF03 <- c(0.445-0.076,0.445+0.076)                                       #Ch10
MgF04 <- c(0.620-0.127,0.620+0.127)                                       #Ch10
MgF05 <- c(0.452-0.086,0.452+0.086)                                       #Ch10
MgF06 <- c(0.693-0.081,0.693+0.081)                                       #Ch10
MgF07 <- c(0.619-0.249,0.619+0.249)                                       #Ch10
MgF08 <- c(0.483-0.183,0.483+0.183)                                       #Ch10
MgF09 <- c(0.546-0.252,0.546+0.252)                                       #Ch10
MgF10 <- c(0.444-0.119,0.444+0.119)                                       #Ch10
MgF11 <- c(0.612-0.131,0.612+0.131)                                       #Ch10
MgF12 <- c(0.465-0.178,0.465+0.178)                                       #Ch10
MgF13 <- c(0.688-0.114,0.688+0.114)                                       #Ch10
MgF14 <- c(0.648-0.175,0.648+0.175)                                       #Ch10
MgF14.5 <- c(0.557-0.172,0.557-0.172) # Average of MgF01 - MgF14
MgF15 <- c(0.30,0.76)                                             #Ch11
MgF16 <- c(0.775-0.055,0.775+0.055) # ≥21 months old Forest-Spring(1 Mar–31 May: 3 months)      #Ch07
MgF17 <- c(0.672-0.012,0.672+0.012) # ≥21 months old Forest-Summer(1 Jun–31 Aug: 3 months)      #Ch07
MgF18 <- c(0.788-0.129,0.788+0.129) # ≥21 months old Forest-Fall-Winter(1 Sep–28 Feb: 6 months) #Ch07
MgF19 <- c(0.881-0.025,0.881+0.025) # ≥21 months old Open-Spring(1 Mar–31 May: 3 months)        #Ch07
MgF20 <- c(0.915-0.085,0.915+0.085) # ≥21 months old Open-Summer(1 Jun–31 Aug: 3 months)        #Ch07
MgF21 <- c(0.931-0.028,0.931+0.028) # ≥21 months old Open-Fall-Winter(1 Sep–28 Feb: 6 months)   #Ch07
MgF21.5 <- c(0.842-0.032,0.842+0.032) # Average of MgF15 - MgF21
MgF22 <- c(0.37-0.15,0.37+0.15)                                           #Ch26
MgF23 <- c(0.58-0.11,0.58+0.11)                                           #Ch23
# MgF24 <- c(0.79) # Incubation Survival                            #Ch25
MgF25 <- c(0.649-0.053,0.649+0.053)                                       #Ch21
# MgF26 <- c(0.53) # Female not Incubate                            #Ch22
# MgF27 <- c(0.47) # Single nest incubation                         #Ch22
MgF28 <- c(0.68-0.02,0.68+0.02)                                           #Ch20
# MgF29 <- c(0.73) # Breeding season survival                       #Ch24


## Reported value for Male Survival

MgM01 <- c(0.604-0.036,0.604+0.036) # ≥21 months old Forest-Spring(1 Mar–31 May: 3 months)      #Ch07
MgM02 <- c(0.960-0.010,0.960+0.010) # ≥21 months old Forest-Summer(1 Jun–31 Aug: 3 months)      #Ch07
MgM03 <- c(0.839-0.065,0.839+0.065) # ≥21 months old Forest-Fall-Winter(1 Sep–28 Feb: 6 months) #Ch07
MgM04 <- c(0.611-0.046,0.611+0.046) # ≥21 months old Open-Spring(1 Mar–31 May: 3 months)        #Ch07
MgM05 <- c(0.848-0.081,0.848+0.081) # ≥21 months old Open-Summer(1 Jun–31 Aug: 3 months)        #Ch07
MgM06 <- c(0.730-0.211,0.730+0.211) # ≥21 months old Open-Fall-Winter(1 Sep–28 Feb: 6 months)   #Ch07
MgM06.5 <- c(0.775-0.068,0.775+0.068)  # Average of MgM01 - MgF06
MgM07 <- c(0.49-0.04,0.49+0.04)                                           #Ch20



