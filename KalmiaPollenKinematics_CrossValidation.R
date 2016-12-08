## Callin Switzer
## 29 Nov 2016
## 8 Dec 2016 Update -- re-digitized pollen and anthers very carefully
## Use cross validation to choose smoothing parameter, using a smoothing spline


## TODO:
# use five -fold cross-validation to decide smoothing parameters or plot noise/resoluation tradeoff
# or simply justify the choice -- by using visual inspection
# compute acceleration and velocity for different values of smoothing parameters
# idea -- show video with smoothed vs. unsmoothed points added -- background subtracted.



# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4", "signal")

ipak(packages)

set.seed(12345)


# read in metadata
dfile <- "/Users/callinswitzer/Dropbox/ExperSummer2015/LaurelsOnly.csv"
metDat <- read.csv(dfile)

metDat <- metDat[metDat$redigitizedFile!= "", ]

# set constants:
fps <- 5000 # frames per second

ii = 3


# read in each .csv file for analysis
# make a list of data frames
digdirect <- "/Users/callinswitzer/Dropbox/ExperSummer2015/CleanKalmiaDigitized/"


newDF <- data.frame()

bestSmthAnth <- numeric()

for(ii in 1:nrow(metDat)){
     
     ddfile <- paste0(digdirect, metDat$redigitizedFile[ii])
     
     dp <- read.csv(ddfile)

     # get anther and pollen locations
     antherPoll <- data.frame(anthx = dp$pt3_cam1_X, anthy= dp$pt3_cam1_Y, 
                              polx = dp$pt4_cam1_X, poly= dp$pt4_cam1_Y)
     # gives only rows where either anth and pollen are complete
     antherPoll = antherPoll[ complete.cases(antherPoll[c('anthx')]) | 
                                   complete.cases(antherPoll[c('polx')]), ]
     
     # if x value starts to right of screen, reverse points,
     # so all x values start on the left part of the screen at 0
     if(lm(antherPoll[,1] ~ I(1:length( antherPoll[,1])))$coefficients[2] < 0 ){
          antherPoll$anthx <- metDat[ii,'vidWidth'] - antherPoll$anthx
          antherPoll$polx <- metDat[ii,'vidWidth'] - antherPoll$polx
     }
     
     
     
     # visualize groups for cv
     
     cvg <- 4
     

     anther <- antherPoll[,c("anthx", "anthy") ]
     anther <- anther[complete.cases(anther), ]
     anther$sampNum <- sample(rep(sample(1:10), length = nrow(anther)))
     anther$time1 <- 1:nrow(anther)
     
     
     newL <- data.frame()
     for(spp in seq(0.001, to = 0.49, by  = 0.01)){
          sumErr <- numeric()
          for(cvg in 1:10){
               xy <- anther
               colnames(xy) <- c("x", "y", "sampeNum", "time1")
               xy$y[anther$sampNum == cvg] <- xy$x[anther$sampNum == cvg]<-  NA
               xy = xy[complete.cases(xy), ]
               #spp <- 0.001
               smy <- smooth.spline(y = xy$y, x = xy$time1, spar = spp)
               smx <- smooth.spline(y = xy$x, x = xy$time1, spar = spp)
               # plot(y = smy$y, x = smx$y, type = 'n')
               # points(xy$x, xy$y, pch = 20) # actual values
               newPredsx = predict(object = smx, x = anther$time1[anther$sampNum == cvg])
               newPredsy = predict(object = smy, x = anther$time1[anther$sampNum == cvg])
               # points(newPredsx$y, y = newPredsy$y, pch = 20, col = 'red') # predicted gaps
               # points(x = anther$anthx[anther$sampNum == cvg], 
               #        y =  anther$anthy[anther$sampNum == cvg],col = 'grey') # actual gaps
               # legend("bottomright", legend = c("raw data", "predicted gaps", "observed gaps"), 
               #        pch = c(20, 20, 1), col = c('black', 'red', 'grey'))
               
               # add up differences -- predicted - observed
               sumdiffs <- data.frame(xObs = anther$anthx[anther$sampNum == cvg],  xPred = newPredsx$y,  
                                      yObs = anther$anthy[anther$sampNum == cvg], yPred = newPredsy$y)
               
               # find distance between points
               sumdiffs$dists = apply(X = sumdiffs, MARGIN = 1, function(x) {
                    return(sqrt((x[1] - x[2])^2 + (x[3] - x[4])^2))
               })
               
               # sum of errors / total number of predicted points
               sumErr[cvg] <- sum(sumdiffs$dists) / nrow(sumdiffs)
          }
          newL <- rbind(newL, c(spp, sum(sumErr) / 10))
     }
     colnames(newL) <- c("smoothingParameter", "sumOfErrors")
     newL
     plot(newL)
     
     bestSmthAnth[ii] <- newL$smoothingParameter[which.min(newL$sumOfErrors)]
}   

hist(bestSmthAnth)
abline(v = mean(bestSmthAnth))
mean(bestSmthAnth) 
median(bestSmthAnth) #,0.291,  maybe go with median, since the dist is skewed

### END OF 10-fold CV for anthers



#### Start CV for Pollen #####
bestSmthPol <- numeric()

for(ii in 1:nrow(metDat)){
     
     ddfile <- paste0(digdirect, metDat$redigitizedFile[ii])
     
     dp <- read.csv(ddfile)
     
     # get anther and pollen locations
     antherPoll <- data.frame(anthx = dp$pt3_cam1_X, anthy= dp$pt3_cam1_Y, 
                              polx = dp$pt4_cam1_X, poly= dp$pt4_cam1_Y)
     # gives only rows where either anth and pollen are complete
     antherPoll = antherPoll[ complete.cases(antherPoll[c('anthx')]) | 
                                   complete.cases(antherPoll[c('polx')]), ]
     
     # if x value starts to right of screen, reverse points,
     # so all x values start on the left part of the screen at 0
     if(lm(antherPoll[,1] ~ I(1:length( antherPoll[,1])))$coefficients[2] < 0 ){
          antherPoll$anthx <- metDat[ii,'vidWidth'] - antherPoll$anthx
          antherPoll$polx <- metDat[ii,'vidWidth'] - antherPoll$polx
     }
     
     
     
     # visualize groups for cv
     pollen <- antherPoll[,c("polx", "poly") ]
     pollen <- pollen[complete.cases(pollen), ]
     pollen$sampNum <- sample(rep(sample(1:10), length = nrow(pollen)))
     pollen$time1 <- 1:nrow(pollen)
     
     
     newL <- data.frame()
     for(spp in seq(0.001, to = 0.49, by  = 0.01)){
          sumErr <- numeric()
          for(cvg in 1:10){
               xy <- pollen
               colnames(xy) <- c("x", "y", "sampeNum", "time1")
               xy$y[pollen$sampNum == cvg] <- xy$x[pollen$sampNum == cvg]<-  NA
               xy = xy[complete.cases(xy), ]
               #spp <- 0.001
               smy <- smooth.spline(y = xy$y, x = xy$time1, spar = spp)
               smx <- smooth.spline(y = xy$x, x = xy$time1, spar = spp)
               # plot(y = smy$y, x = smx$y, type = 'n')
               # points(xy$x, xy$y, pch = 20) # actual values
               newPredsx = predict(object = smx, x = pollen$time1[pollen$sampNum == cvg])
               newPredsy = predict(object = smy, x = pollen$time1[pollen$sampNum == cvg])
               # points(newPredsx$y, y = newPredsy$y, pch = 20, col = 'red') # predicted gaps
               # points(x = pollen$polx[pollen$sampNum == cvg],
               #        y =  pollen$poly[pollen$sampNum == cvg],col = 'grey') # actual gaps
               # legend("bottomright", legend = c("raw data", "predicted gaps", "observed gaps"), 
               #        pch = c(20, 20, 1), col = c('black', 'red', 'grey'))
               
               # add up differences -- predicted - observed
               sumdiffs <- data.frame(xObs = pollen$polx[pollen$sampNum == cvg],  xPred = newPredsx$y,  
                                      yObs = pollen$poly[pollen$sampNum == cvg], yPred = newPredsy$y)
               
               # find distance between points
               sumdiffs$dists = apply(X = sumdiffs, MARGIN = 1, function(x) {
                    return(sqrt((x[1] - x[2])^2 + (x[3] - x[4])^2))
               })
               
               # sum of errors / total number of predicted points
               sumErr[cvg] <- sum(sumdiffs$dists) / nrow(sumdiffs)
          }
          newL <- rbind(newL, c(spp, sum(sumErr) / 10))
     }
     colnames(newL) <- c("smoothingParameter", "sumOfErrors")
     newL
     plot(newL)
     
     bestSmthPol[ii] <- newL$smoothingParameter[which.min(newL$sumOfErrors)]
}   

hist(bestSmthPol)
abline(v = mean(bestSmthPol)) # 
median(bestSmthPol) # 0.286, maybe go with median, since the dist is skewed
mean(bestSmthPol)#  0.25725


# Decision: smoothing parameter is between 0.286 and 0.296 -- go with 0.29

###########


# Visualize raw and smoothed data to see how different it is.

spp = 0.29
for(ii in 1:nrow(metDat)){
     
     ddfile <- paste0(digdirect, metDat$redigitizedFile[ii])
     
     dp <- read.csv(ddfile)
     
     # get anther and pollen locations
     antherPoll <- data.frame(anthx = dp$pt3_cam1_X, anthy= dp$pt3_cam1_Y, 
                              polx = dp$pt4_cam1_X, poly= dp$pt4_cam1_Y)
     # gives only rows where either anth and pollen are complete
     antherPoll = antherPoll[ complete.cases(antherPoll[c('anthx')]) | 
                                   complete.cases(antherPoll[c('polx')]), ]
     
     # if x value starts to right of screen, reverse points,
     # so all x values start on the left part of the screen at 0
     if(lm(antherPoll[,1] ~ I(1:length( antherPoll[,1])))$coefficients[2] < 0 ){
          antherPoll$anthx <- metDat[ii,'vidWidth'] - antherPoll$anthx
          antherPoll$polx <- metDat[ii,'vidWidth'] - antherPoll$polx
     }
     
     
     
     # visualize groups for cv
     pollen <- antherPoll[,c("polx", "poly") ]
     pollen <- pollen[complete.cases(pollen), ]
     pollen$sampNum <- sample(rep(sample(1:10), length = nrow(pollen)))
     pollen$time1 <- 1:nrow(pollen)
     
     
     
     xy <- pollen
     colnames(xy) <- c("x", "y", "sampeNum", "time1")
     xy = xy[complete.cases(xy), ]
     smy <- smooth.spline(y = xy$y, x = xy$time1, spar = spp)
     smx <- smooth.spline(y = xy$x, x = xy$time1, spar = spp)
     plot(y = smy$y, x = smx$y, pch = 20, col = rgb(0,0,0,0.5), asp = 1)
     points(xy$x, xy$y, pch = 20, col = rgb(1,0,0, 0.5)) # actual values
     
     Sys.sleep(0.5)
}   


spp = 0.29
for(ii in 1:nrow(metDat)){
     
     ddfile <- paste0(digdirect, metDat$redigitizedFile[ii])
     
     dp <- read.csv(ddfile)
     
     # get anther and pollen locations
     antherPoll <- data.frame(anthx = dp$pt3_cam1_X, anthy= dp$pt3_cam1_Y, 
                              polx = dp$pt4_cam1_X, poly= dp$pt4_cam1_Y)
     # gives only rows where either anth and pollen are complete
     antherPoll = antherPoll[ complete.cases(antherPoll[c('anthx')]) | 
                                   complete.cases(antherPoll[c('polx')]), ]
     
     # if x value starts to right of screen, reverse points,
     # so all x values start on the left part of the screen at 0
     if(lm(antherPoll[,1] ~ I(1:length( antherPoll[,1])))$coefficients[2] < 0 ){
          antherPoll$anthx <- metDat[ii,'vidWidth'] - antherPoll$anthx
          antherPoll$polx <- metDat[ii,'vidWidth'] - antherPoll$polx
     }
     
     
     
     # visualize groups for cv
     anther <- antherPoll[,c("anthx", "anthy") ]
     anther <- anther[complete.cases(anther), ]
     anther$time1 <- 1:nrow(anther)
     
     
     
     xy <- anther
     colnames(xy) <- c("x", "y", "time1")
     xy = xy[complete.cases(xy), ]
     smy <- smooth.spline(y = xy$y, x = xy$time1, spar = spp)
     smx <- smooth.spline(y = xy$x, x = xy$time1, spar = spp)
     plot(y = smy$y, x = smx$y, pch = 20, col = rgb(0,0,0,0.5), asp = 1)
     points(xy$x, xy$y, pch = 20, col = rgb(1,0,0, 0.5)) # actual values
     legend("bottomright", pch = c(20,20), col = c('black', 'red'), legend = c("predicted", "observed"))
     
     Sys.sleep(0.5)
}   






