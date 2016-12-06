## Laurel Digitizations
# Callin Switzer
# June 21, 2015
# Updated 16 Nov 2016 - redone to make code better
# Calculate velocity, position, and acceleration, of anther and pollen release 
# from mountain laurels from Arboretum.

## TODO:
# smooth, and impute points to get better estimates
# try different smoothing methods -- wavelet analysis?, butterworth filter?



# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4")

ipak(packages); rm(packages); rm(ipak)


# read in metadata
dfile <- "/Users/callinswitzer/Dropbox/ExperSummer2015/LaurelsOnly.csv"
metDat <- read.csv(dfile)

metDat <- metDat[metDat$digitizedFile!= "", ]

# set constants:
fps <- 5000 # frames per second

# read in each .csv file for analysis
# make a list of data frames
digdirect <- "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/"



# use ii = 2, to make an example movie with digitized, smoothed pollen and anther
# 20150608_135726.avi is the movie, where ii = 2
# called movie digPollVid_20150608_135726.mp4
newL = list()
for(ii in 1:length(metDat$digitizedFile)){
     ddfile <- paste0(digdirect, metDat$digitizedFile[ii])
     
     dp <- read.csv(ddfile)
     
     # calibrate locations, based on digitized pin or other object
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     
     # get the number of pixels in the calibration
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X - pin$dp.pt2_cam1_X)^2 +
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / 
          metDat$CalSizeMM[ii] # to get to mm
     
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
     
     # #play video
     # plot(c(antherPoll$polx,antherPoll$anthx), 
     #      c(antherPoll$poly, antherPoll$anthy), 
     #      type = 'n', asp = 1)
     #      
     # # slow method, but can make a movie, if needed
     # for(jj in 1:nrow(antherPoll)){
     #      points(antherPoll$polx[jj], antherPoll$poly[jj], pch = 20)
     #      points(antherPoll$anthx[jj], antherPoll$anthy[jj], pch = 20, col = 'red')
     #     # Sys.sleep(0.1)
     #      
     # }
     
     
     # cbind data frame, to add smoothed columns
     antherPoll <- data.frame(cbind(antherPoll, antherPoll))
     
     
     # smooth
     smspar <- 0.5 # smoothing parameter
     foo = sapply(X = c("anthx.1", "anthy.1", "polx.1", "poly.1"), FUN = function(x){
          sm1 <- smooth.spline(na.omit(antherPoll[, x]), spar = smspar) 
          antherPoll[, x][complete.cases(antherPoll[, x])] <<- sm1$y
     })

     # # check smooth
     # plot(c(antherPoll$polx,antherPoll$anthx),
     #      c(antherPoll$poly, antherPoll$anthy),
     #      type = 'n', asp = 1)
     # 
     # points(antherPoll$polx, antherPoll$poly, pch = 20)
     # points(antherPoll$anthx, antherPoll$anthy, pch = 20, col = 'red')
     # lines(antherPoll$polx.1, antherPoll$poly.1, col = 'grey', type = 'b', pch = 20, cex = 0.5)
     # lines(antherPoll$anthx.1, antherPoll$anthy.1,col = 'blue', type = 'b', pch = 20, cex = 0.5)
     
     # write an example .csv file to make into a video --
     # write.csv("~/Desktop/DigKalmPts.csv", x = antherPoll, row.names = FALSE)
     
     # add time to data frame
     antherPoll$tme = 0: (nrow(antherPoll) - 1) / fps # time
     
     # add columns with absolute position into dataframe (calculated from smoothed data)
     # calculate position from starting point, not from minimum point
     bar = sapply(X = c("anthx.1", "anthy.1", "polx.1", "poly.1"), FUN = function(x){
          newName = paste0(x, ".abs")
          tmp <- antherPoll[,x] / PixInPin / 1000
          antherPoll[,newName] <<- tmp - na.omit(tmp)[1]
          #antherPoll[,newName] <<- tmp - min(na.omit(tmp))
     })
     
     # # check position
     # plot(c(antherPoll$polx.1.abs,antherPoll$anthx.1.abs),
     #      c(antherPoll$poly.1.abs, antherPoll$anthy.1.abs), asp = 1)
     
     
     # add columns to show velocity, based on smoothed, absolute position
     # velocity is in m/s
     bat = sapply(X = c("anthx.1.abs", "anthy.1.abs", "polx.1.abs", "poly.1.abs"), 
                  FUN = function(x){
          newName = paste0(x, ".vel")
          tmp <-  c(NaN, diff(antherPoll[,x])) * fps # add a NaN to beginning of data
          antherPoll[,newName] <<- tmp 
     })
     
     
     # add columns to show velocity magnitude
     # adds X and Y velocities, using Pythagorean Theorem
     antherPoll$anthVelMag = sqrt(antherPoll$anthx.1.abs.vel^2 + antherPoll$anthy.1.abs.vel^2)
     antherPoll$polVelMag = sqrt(antherPoll$polx.1.abs.vel^2 + antherPoll$poly.1.abs.vel^2)
     
     # calculate acceleration
     bam = sapply(X = c("anthx.1.abs.vel", "anthy.1.abs.vel", "polx.1.abs.vel", "poly.1.abs.vel"), 
                   FUN = function(x){
                        newName = paste0(x, ".acc")
                        # add a NaN to beginning of data
                        tmp <-  c(NaN, diff(antherPoll[,x])) * fps 
                        antherPoll[,newName] <<- tmp 
                   })
     
     #position
     plot(antherPoll$anthx.1.abs, y = antherPoll$anthy.1.abs, type = 'b', ylim = c(0, 0.06), asp = 1)
     lines(antherPoll$polx.1.abs, antherPoll$poly.1.abs, type = 'b', xlim = c(0, 50), col = 'red')
     
     par(mfrow = c(3,1), mai = c(0,.5,0,0))
     plot(antherPoll$anthx.1.abs, type = 'l', xlim = c(0, 50), ylim = c(-0.004, 0.018), ylab = "anther position")
     lines(antherPoll$anthy.1.abs, type = 'l', xlim = c(0, 50), col = 'blue')
     lines(antherPoll$distanth, type = 'l', col = 'red')
     legend("topright", legend = c("x", "y", "total"), lty = 1, col = c("black", "blue", "red"))
     abline(h = 0, lty = 2)
     abline(v = c(15, 14, 13), lty = 2)
     
     
     # velocity
     plot(antherPoll$anthx.1.abs.vel, type = 'l', xlim = c(0,50), ylim = c(-3.1, 4.1), ylab = "anther speed")
     lines(antherPoll$anthy.1.abs.vel, type = 'l', col = 'blue')
     lines(antherPoll$anthVelMag, type = 'l', col = 'red')
     abline(h = 0, lty = 2)
     abline(v = c(15, 14, 13), lty = 2)
     
     
     plot(antherPoll$anthx.1.abs.vel.acc, type = 'l', xlim = c(0, 50), ylab = "anther acceleration", ylim =c(-3000, 3000))
     lines(antherPoll$anthy.1.abs.vel.acc, type = 'l', xlim = c(0, 50), col = 'blue')
     lines(antherPoll$anthAcc.1)
     a1 <- c(NA, diff(antherPoll$anthVelMag))*fps
     lines(a1, col = 'red')
     abline(h = 0, lty = 2)
     abline(v = c(15, 14, 13), lty = 2)
     

     
     # add columns to show acceleration -- should be in m/s/s
     # adds X and Y accelerations, using Pythagorean Theorem
    
    signedSqAcc = antherPoll$anthx.1.abs.vel.acc^2  + 
                   antherPoll$anthy.1.abs.vel.acc^2  
     
     antherPoll$anthAcc.1 <- sqrt(abs(signedSqAcc)) * sign(signedSqAcc)
     
     signedSqAcc_pol = antherPoll$polx.1.abs.vel.acc^2 * 
                         sign(antherPoll$polx.1.abs.vel.acc) +
                              antherPoll$poly.1.abs.vel.acc * 
                         sign(antherPoll$poly.1.abs.vel.acc)
     antherPoll$polAcc.1 <- sqrt(abs(signedSqAcc_pol)) * sign(signedSqAcc_pol)
     
     
     # calculate displacement (total distance, in meters)
     dd1 <- sqrt(diff(antherPoll$anthx.1.abs)^2 + 
                      diff(antherPoll$anthy.1.abs)^2)
     dd1[is.na(dd1)] <- 0
     antherPoll$distanth <- c(0, cumsum(dd1))
     antherPoll$distanth[is.na(antherPoll$anthy.1.abs)] <- NaN
     
     dd2 <- sqrt(diff(antherPoll$polx.1.abs)^2 + 
                      diff(antherPoll$poly.1.abs)^2)
     dd2[is.na(dd2)] <- 0
     antherPoll$distpoll <- c(0, cumsum(dd2))
     antherPoll$distpoll[is.na(antherPoll$poly.1.abs)] <- NaN
     
     
     # #check
     # plot(antherPoll$polx.1.abs, antherPoll$poly.1.abs, asp = 1)
     # 
     # with(antherPoll, {
     #      plot(y = distpoll,x =  tme)
     #      points( y=poly.1.abs, x = tme, pch = 20, col = rgb(0,1,1,0.5))
     # })
     # 
     # plot(antherPoll$anthx.1.abs, antherPoll$anthy.1.abs, asp = 1)
     # 
     # with(antherPoll, {
     #      plot(y = antherPoll$distanth, x =  antherPoll$tme, col = rgb(1,0,0,1))
     #      points( y=anthy.1.abs, x = tme, pch = 20, col = rgb(0,1,0,0.5))
     #      points( y=anthx.1.abs, x = tme, pch = 20, col = rgb(0,1,0,0.5))
     # })
     
     # # smooth
     # smspar <- 0.5 # smoothing parameter
     # foob = sapply(X = c("anthx.1.abs", "anthy.1.abs", 
     #                     "polx.1.abs", "poly.1.abs"), FUN = function(x){
     #      
     #      # fit smooth
     #      sm1 <- smooth.spline(y = na.omit(antherPoll[, x]), 
     #                           x = antherPoll$tme[!is.na(antherPoll[, x])], 
     #                           spar = smspar) 
     #      
     #      # predict new data via interpolation
     #      sm2 <- predict(sm1, x = seq(min(antherPoll$tme[!is.na(antherPoll[, x])]), 
     #                                  max(antherPoll$tme[!is.na(antherPoll[, x])]), 
     #                                  length.out = 1000))
     #      
     #      newName <- paste0(strsplit(x, "\\.")[[1]][1], "_disp")
     #      newT <- paste0(newName, "_time")
     #      
     #      return(data.frame(sm2))
     #      
     # })
     # 
     # foob[['x']]
     # 
     
     
     # add column to show total distance (from minimum) 
     # in meters
     # this is not the same as total distance traveled
     # antherPoll$polDispMag = sqrt(antherPoll$polx.1.abs^2 +
     #                                antherPoll$poly.1.abs^2)
     # antherPoll$anthDispMag = sqrt(antherPoll$anthx.1.abs^2 + 
     #                                   antherPoll$anthy.1.abs^2)
     # 
     # plot(antherPoll$polDispMag)
     # points(antherPoll$distpoll,pch = 20,  col = rgb(0,1,0, 0.5))
     # plot(antherPoll$anthDispMag, ylim = c(0,0.02))
     # points(antherPoll$distanth,pch = 20,  col = rgb(0,1,0, 0.5))
     
     # check
     # plot(y = antherPoll$anthAccMag,  x= antherPoll$tme, xlim = c(0, 0.02), type = 'l')
     
     # Sys.sleep(0.1)
     print(ii)
     
     smList = list(file = ddfile, 
                   pixelConversion = PixInPin, 
                   maxAnthVel = max(antherPoll$anthVelMag, na.rm = TRUE), 
                   maxPollVel = max(antherPoll$polVelMag, na.rm = TRUE), 
                   maxAnthAcc = max(antherPoll$anthAcc.1, na.rm = TRUE), 
                   maxPollAcc = max(antherPoll$polAcc.1, na.rm = TRUE), 
                   anthPollDF = antherPoll)
     newL[[ii]] = smList
}
#  ii = 21, 25, 28: pollen digitization may need smoothing
# check: 20150608_143846 (ii = 5)

newL[[1]]$file
with(newL[[1]], c(maxAnthVel, maxPollVel))

ii = 1
plot(newL[[ii]]$anthPollDF$distanth)

# save dataframe
saveDir <- "/Users/callinswitzer/Dropbox/ExperSummer2015/Kalmia2015FiguresAndData/"
save(newL, file = paste0(saveDir, "Kalm2015DigitizedDataset.rda"))

# how to load in dataset, if you need to do it in the future
# load(paste0(saveDir, "Kalm2015DigitizedDataset.rda"))


##################################################################
##
## Make Velocity, Acceleration, and Displacement Curves
##
##################################################################

#[REFREF]
### HERE: Figure out displacement -- want to align with steep slopes at 0
## Why isn't velocity max when slope is greatest
# calculate velocity in two different ways and see if they're the same!

ii = 1
plot(y = newL[[ii]]$anthPollDF$distanth, x = newL[[ii]]$anthPollDF$tme, 
     xlab = 'time (s)', ylab = 'displacement from initial position (m)', 
     xlim = c(-0.012, 0.012), ylim = c(0,0.08), type = 'n')

# align sections, based on max velocity
for(ii in 1:length(newL)){
     tmp2 <- newL[[ii]]$anthPollDF
     timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$anthVelMag)]
     
     lines(y = newL[[ii]]$anthPollDF$distanth, x = timeCentered, 
           col = 'red')
}
abline(v = 0)


# velocity
plot(y = newL[[ii]]$anthPollDF$anthVelMag, x = newL[[ii]]$anthPollDF$tme, 
     xlab = 'time (s)', ylab = 'velocity (m/s)', type = 'n', 
     xlim = c(-0.01, 0.02), ylim = c(0, 6))

# align sections, based on max velocity
for(ii in 1:length(newL)){
     tmp2 <- newL[[ii]]$anthPollDF
     timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$anthVelMag)]
     
     lines(y = newL[[ii]]$anthPollDF$anthVelMag, x = timeCentered, 
           col = 'red')
}
abline(v = 0)


# calculate velocity in a different way, and see if it's the same
# velocity
plot(y = newL[[ii]]$anthPollDF$anthVelMag, x = newL[[ii]]$anthPollDF$tme, 
     xlab = 'time (s)', ylab = 'velocity (m/s)', type = 'n', 
     xlim = c(-0.01, 0.02), ylim = c(0, 6))

# align sections, based on max velocity
# should be called speed, since it doesn't take direction into account
for(ii in 1:length(newL)){
     tmp2 <- newL[[ii]]$anthPollDF
     timeCentered <- tmp2$tme- tmp2$tme[which.max(tmp2$anthVelMag)]
     
     vel <- c(NA, diff(tmp2$distanth)*5000)
     
     lines(y = vel, x = timeCentered, 
           col = 'green')
}
abline(v = 0) # same


# accel
plot(y = newL[[ii]]$anthPollDF$anthAcc.1, x = newL[[ii]]$anthPollDF$tme, 
     xlab = 'time (s)', ylab = 'acceleration (m/s/s)', type = 'n', 
     xlim = c(-0.01, 0.02), ylim = c(-6000, 6000))

# align sections, based on max velocity
for(ii in 1:length(newL)){
     tmp2 <- newL[[ii]]$anthPollDF
     timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$anthAcc.1)]
     
     lines(y = newL[[ii]]$anthPollDF$anthAcc.1, x = timeCentered, 
           col = 'red')
}
abline(v = 0)


# calculate acc in a different way, and see if it's the same
# acceleration
plot(y = newL[[ii]]$anthPollDF$anthAcc.1, x = newL[[ii]]$anthPollDF$tme, 
     xlab = 'time (s)', ylab = 'acceleration (m/s/s)', type = 'n', 
     xlim = c(-0.01, 0.02), ylim = c(-3000, 3000))

# align sections, based on max velocity
for(ii in 1:length(newL)){
     tmp2 <- newL[[ii]]$anthPollDF
     timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$anthAcc.1)]
     
     v1 <- diff(tmp2$distanth)*5000
     acc <- c(NA, NA, (diff(v1)*5000))
     
     lines(y = acc, x = timeCentered, 
           col = 'pink')
     abline(h = 0)
     
     lines(y = tmp2$anthAcc.1, x = timeCentered, 
           col = 'red')
}
abline(v = 0) # not the same

###############################################
### LOOK AT DIST, VEL, ACC

par(mfrow = c(3,1))
ii = 1
tmp2 <- newL[[ii]]$anthPollDF


#plot distance traveled for anther
timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$polVelMag)]
plot(y = tmp2$distpoll, x = timeCentered, 
      col = 'red')

# velocity
par(mfrow= c(2,1))
timeCentered <- tmp2$tme- tmp2$tme[which.max(tmp2$anthVelMag)]
vel <- c(NA, diff(tmp2$distanth)*5000)
plot(y = vel, x = timeCentered, 
      col = 'green', type = 'l', xlim = c(0, 0.01))


# Acceleration
timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$anthVelMag)]

v1 <- diff(tmp2$distanth)*5000
acc <- c(NA, NA, (diff(v1)*5000))

plot(y = acc, x = timeCentered, 
      col = 'pink', type = 'l', xlim = c(0, 0.01), ylim = c(-3000, 3000))
lines(y = tmp2$anthAcc.1, x = timeCentered, 
      col = 'blue', pch = 20)
abline(h = 0)





timeCentered <- tmp2$tme - tmp2$tme[which.max(tmp2$anthAcc.1)]

v1 <- diff(tmp2$distanth)*5000
acc <- c(NA, NA, (diff(v1)*5000))

lines(y = acc, x = timeCentered, 
      col = 'pink')
abline(h = 0)

lines(y = newL[[ii]]$anthPollDF$anthAcc.1, x = timeCentered, 
      col = 'red')





##### OLD BELOW THIS LINE ######################



#Anther Video
ii = 31
anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
anther <- anther[ complete.cases(anther), ]
anther$x <- anther$x - min(anther$x)
anther$y <- anther$y - min(anther$y)

anther <- anther / PixInPin / 1000 # in meters
par(pty="s")

plot(anther$x, anther$y, ylab = "y position (meters)", xlab = "x position (meters)", 
                 xlim = c(0, 0.01), ylim = c(0, 0.01), pch = ".", col = rgb(0,0,0, 0.1))
for(j in 1:length(anther$x)){
     points(anther$x[j], anther$y[j], pch = 20, col= rgb(0,0,0, 0.2), type = 'p')
     Sys.sleep(0.2)
}

smspar <- 0.5
xsmth <- smooth.spline(anther$x, spar = smspar) 
ysmth <- smooth.spline(anther$y, spar = smspar)

lines(xsmth$y, ysmth$y, type = "l", lwd = 3, col = rgb(0, 0, 0, 0.5)) 


tme <- 1:length(anther$x) / fps # time

# Smooth

smspar <- 0.5
xsmth <- smooth.spline(anther$x, spar = smspar) 
ysmth <- smooth.spline(anther$y, spar = smspar)
vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))
velMag <- with(vel, sqrt(xvel^2 + yvel^2))*5000

indMax <- which(velMag == max(velMag))

ttme <- -indMax:(length(velMag) - indMax) / 5000
par(pty="m")

plot(y = abs(vel$xvel), x = ttme[-length(ttme)], pch = ".", xlab = "time (s)", ylab = "x velocity (m/s)", type = "l")

plot(y = abs(vel$yvel), x = ttme[-length(ttme)], pch = ".", xlab = "time (s)", ylab = "y velocity (m/s)", type = "l")

plot(y = velMag, x = ttme[-length(ttme)], pch = ".", xlab = "time (s)", ylab = "total velocity (m/s)", type = "n")
lines(y = velMag, x = ttme[-length(tme)], col = rgb(0,0,0,.5), lwd = 2)


# velocity
figDir = "/Users/callinswitzer/Dropbox/ExperSummer2016/Kalmia/KalmiaFigures/"

par(pty="m")
pdf(file = paste0(figDir, "velocityTrace.pdf"), width = 5, height = 4)

velMaxAnth <- numeric(length(metDat$digitizedFile))
for(ii in 1:length(metDat$digitizedFile)){
     ddfile <- paste("/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", metDat$digitizedFile[ii], sep = "")
     
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] # to get to mm
     
     #Anther
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     
     anther <- anther / PixInPin / 1000 # in meters
     
#      if(ii == 1) plot(anther$x, anther$y, ylab = "y position (meters)", xlab = "X position (meters)", 
#                       xlim = c(0, 0.015), ylim = c(0, 0.015), pch = ".")
#      points(anther$x, anther$y, pch = ".")
     tme <- 1:length(anther$x) / fps # time
     
     
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
#      lines(xsmth$y, ysmth$y, type = "l", col = rgb(0, 0, 0, 0.5))

     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))


     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))*5000

     indMax <- which(velMag == max(velMag))

ttme <- -indMax:(length(velMag) - indMax) / 5000


     
     if (ii == 1) plot(y = velMag, x = ttme[-length(ttme)], pch = ".", xlab = "time (s)", ylab = "velocity (m/s)", xlim = c(-0.012, 0.012), ylim = c(0,6), type = "n")

     lines(y = velMag, x = ttme[-length(tme)], col = rgb(0,0,0,0.5))
     
     Sys.sleep(0.2)
     
     velMaxAnth[ii] <- velMag[indMax]

     
}

dev.off()



# Pollen Velocity
par(pty="m")
pdf(file = paste0(figDir, "velocityTracePollen.pdf"), width = 5, height = 4)

velMaxPol <- numeric(length(metDat$digitizedFile))
for(ii in 1:length(metDat$digitizedFile)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #Pollen
     pollen <- data.frame(x = dp$pt4_cam1_X, y= dp$pt4_cam1_Y)
     pollen <- pollen[ complete.cases(pollen), ]
     pollen$x <- pollen$x - min(pollen$x)
     pollen$y <- pollen$y - min(pollen$y)
     
     pollen <- pollen / PixInPin / 1000 # in meters

     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(pollen$x, spar = smspar) 
     ysmth <- smooth.spline(pollen$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))*5000
     indMax <- which(velMag == max(velMag))
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if (ii == 1){
          plot(y = velMag, x = ttme[-length(ttme)], pch = ".", 
               xlab = "time (s)", ylab = "velocity (m/s)", 
               xlim = c(-0.012, 0.1), ylim = c(0,5), type = "n")  
     } 
     
     lines(y = velMag, x = ttme[-length(ttme)], col = rgb(1,0,.3,0.5))
     velMaxPol[ii] <- velMag[indMax]    
}

dev.off()


# LM
regMod <- lm(velMaxAnth~1)
confint(regMod) # too wide!

#LMER

modVelMax <- lmer(formula = velMaxAnth ~ (1|plant/FlowerNumber), data = metDat)
summary(modVelMax)
confint(modVelMax)
plot(modVelMax)


# BS CI for the maximum (without dependencies)
bsFunc <- function(dat){
     mean(sample(velMaxAnth, replace = TRUE))
}

foo <- replicate(n = 10000, bsFunc(metDat))

hist(foo)
quantile(foo, probs = c(0.025, 0.975))




# acceleration
#accel <- 





#pollen
pollen <- data.frame(x = dp$pt4_cam1_X, y = dp$pt4_cam1_Y)
pollen <- pollen[ complete.cases(pollen), ]
pollen <- pollen / PixInPin / 100 # in meters
plot(y =  pollen$y, x = pollen$x, ylim = c(0,.10), xlim = c(0, .10))

plot(pollen$x, pollen$y, ylab = "meters", xlab = )
tme <- 1:length(pollen$x) / fps

smspar <- .5

xsmth <- smooth.spline(pollen$x, spar = smspar) 

ysmth <- smooth.spline(pollen$y, spar = smspar)

lines(xsmth$y, ysmth$y, type = "l", lwd = 4, col = "red")


# velocity
vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))

velMag <- with(vel, sqrt(xvel^2 + yvel^2))*5000

plot(y = velMag, x = tme[-length(tme)])

####################################################################################
## Acceleration

par(pty="m")
par(mfrow = c(3,1))
for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
                     metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 

     #ANTHER
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     anther <- anther / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
#      acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000

     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if(ii == 1){
          plot(y = ysmth$y, x= ttme, type = "l", 
               col = rgb(0,0,0, 0.5), ylim = c(0, .015), xlim = c(-0.008, 0.012),  
               xlab = "time (s)", ylab = "y position (m)")
     } 
     else lines(y = ysmth$y, x= ttme, col = rgb(0,0,0, 0.5))
#      plot(vel$yvel, type = "l")
#      plot(acc$yacc, type = "l")
}

for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #ANTHER
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     anther <- anther / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     #      acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if(ii == 1){
          plot(y = vel$yvel, x= ttme[-length(ttme)], type = "l", 
               col = rgb(1,0,0.5, 0.3), xlim = c(-0.008, 0.012), ylim = c(-5, 6),
               xlab = "time (s)", ylab = "y velocity (m/s)")
     } 
     else lines(y = vel$yvel, x= ttme[-length(ttme)], col = rgb(1,0,0.3, 0.5))
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}


for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #ANTHER
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     anther <- anther / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if(ii == 1){
          plot(y = acc$yacc, x= ttme[1:(length(ttme)-2)], type = "l", 
               col = rgb(0,0.3,1, 0.5),  xlim = c(-0.008, 0.012), ylim = c(-5000, 5000),
               xlab = "time (s)", ylab = "y acceleration (m/s/s)")
     } 
     else lines(y = acc$yacc, x= ttme[1:(length(ttme)-2)], col = rgb(0,0.3,1, 0.5))
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}



#########################################################
### X
par(pty="m")
par(mfrow = c(3,1))
for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #ANTHER
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     # center anther, and make sure all anthers are traveling in 
     # the same direction
     anther$x <- scale(anther$x, center = TRUE, scale = FALSE)
     if(lm(anther$x ~ I(1:length(anther$x)))$coefficients[2] > 0){
          anther$x <- -anther$x
     }
     
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     anther <- anther / PixInPin / 1000 # in meters
     
#     if(anther$x[1] > 0.005) anther$x <- rev(anther$x)
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     #      acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if(ii == 1){
          plot(y = xsmth$y, x= ttme, type = "l", 
               col = rgb(0,0,0, 0.5),  xlim = c(-0.008, 0.012),  ylim = c(0, 0.02), 
               xlab = "time (s)", ylab = "x position (m)")
     } 
     else lines(y = xsmth$y, x= ttme, col = rgb(0,0,0, 0.5))
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}

for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #ANTHER
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     # center anther, and make sure all anthers are traveling in 
     # the same direction
     anther$x <- scale(anther$x, center = TRUE, scale = FALSE)
     if(lm(anther$x ~ I(1:length(anther$x)))$coefficients[2] > 0){
          anther$x <- -anther$x
     }
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     anther <- anther / PixInPin / 1000 # in meters
     
     if(anther$x[1] > 0.005) anther$x <- rev(anther$x)
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     #      acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if(ii == 1){
          plot(y = vel$xvel, x= ttme[-length(ttme)], type = "l", 
               col = rgb(1,0,0.5, 0.3), xlim = c(-0.008, 0.012), ylim = c(-5, 6),
               xlab = "time (s)", ylab = "x velocity (m/s)")
     } 
     else lines(y = vel$xvel, x= ttme[-length(ttme)], col = rgb(1,0,0.3, 0.5))
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}


for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #ANTHER
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     # center anther, and make sure all anthers are traveling in 
     # the same direction
     anther$x <- scale(anther$x, center = TRUE, scale = FALSE)
     if(lm(anther$x ~ I(1:length(anther$x)))$coefficients[2] > 0){
          anther$x <- -anther$x
     }
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     anther <- anther / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     if(ii == 1){
          plot(y = acc$xacc, x= ttme[1:(length(ttme)-2)], type = "l", 
               col = rgb(0,0.3,1, 0.5),  xlim = c(-0.008, 0.012), ylim = c(-5000, 5000),
               xlab = "time (s)", ylab = "x acceleration (m/s/s)")
     } 
     else lines(y = acc$xacc, x= ttme[1:(length(ttme)-2)], col = rgb(0,0.3,1, 0.5))
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}


##############################################################################
## total accel
par(pty="m")
pdf(file = "accelTrace.pdf", width = 5, height = 4)

accMaxAnth <- numeric(nrow(metDat))

for(ii in 1:length(metDat$digitizedFile)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #Anther
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther <- anther[ complete.cases(anther), ]
     # center anther, and make sure all anthers are traveling in 
     # the same direction
     anther$x <- scale(anther$x, center = TRUE, scale = FALSE)
     if(lm(anther$x ~ I(1:length(anther$x)))$coefficients[2] > 0){
          anther$x <- -anther$x
     }
     anther$x <- anther$x - min(anther$x)
     anther$y <- anther$y - min(anther$y)
     
     anther <- anther / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))*5000
     acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))
     
     accMag <- with(acc, sqrt(xacc^2 + yacc^2))*5000^2
     #plot(accMag, type = "l")
     indMax <- which(velMag == max(velMag))
     indMaxAcc <- which(accMag == max(accMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     tmeAcc <- -indMaxAcc:(length(accMag) - indMaxAcc) / 5000
     
     #plot(y = accMag, x = tmeAcc[-length(tmeAcc)], type = 'l')

     
     
     
     if (ii == 1) plot(y = accMag, x = tmeAcc[-length(tmeAcc)], pch = ".", xlab = "time (s)", ylab = "acceleration (m/s/s)", xlim = c(-0.012, 0.012), ylim = c(0, 6000), type = "n")
     
     lines(y = accMag, x = tmeAcc[-length(tmeAcc)], col = rgb(0,0,0,0.5))
     
     Sys.sleep(0.2)
     
     accMaxAnth[ii] <- accMag[indMaxAcc]
     
     
}

dev.off()

# LMER for acceleration
modAccMax <- lmer(formula = log(accMaxAnth) ~ (1|plant/FlowerNumber), data = metDat)
plot(modAccMax)
exp(confint(modAccMax))

summary(modAccMax)
plot(modAccMax)
bootF <- function(.){
     c(beta = fixef(.))
}

boo1 <- bootMer(modAccMax, bootF, nsim = 1000)
require(boot)
hist(boo1$t)
# CI
(bCI.1 <- boot.ci(boo1, index=1, type=c("norm", "basic", "perc")))# beta
exp(bCI.1$normal)


####################################################################################
# Pollen total Y Distance
MaxHeight <- numeric()
for(ii in 1:length(metDat$digitizedFile)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
#      #Anther
#      anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
#      anther <- anther[ complete.cases(anther), ]
#      anther$x <- anther$x - min(anther$x)
#      anther$y <- anther$y - min(anther$y)
     
#      anther <- anther / PixInPin / 1000 # in meters

     #Pollen
     pollen <- data.frame(x = dp$pt4_cam1_X, y= dp$pt4_cam1_Y)
     pollen <- pollen[ complete.cases(pollen), ]
     pollen$x <- pollen$x - min(pollen$x)
     pollen$y <- pollen$y - min(pollen$y)
     pollen <- pollen / PixInPin / 1000 # in meters

     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(pollen$x, spar = smspar) 
     ysmth <- smooth.spline(pollen$y, spar = smspar)
     
     #plot(ysmth)

     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))*5000
     acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))
     
     accMag <- with(acc, sqrt(xacc^2 + yacc^2))*5000^2
     #plot(accMag, type = "l")
     indMax <- which(velMag == max(velMag))
     indMaxAcc <- which(accMag == max(accMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     tmeAcc <- -indMaxAcc:(length(accMag) - indMaxAcc) / 5000
     
     #plot(y = accMag, x = tmeAcc[-length(tmeAcc)], type = 'l')
     
     
     
     
     if (ii == 1) plot(y = pollen$y, x = ttme, pch = ".", xlab = "time (s)", ylab = "pollen height (m)", xlim = c(0, 0.1), ylim = c(0, 0.1), type = "n")
     
     lines(y = pollen$y, x = ttme, col = rgb(0,0,0,0.5))
     
     Sys.sleep(0.2)
     
     MaxHeight[ii] <- max(pollen$y)
     
     
}


# LMER for max height
modAccMax <- lmer(formula = log(MaxHeight) ~ (1|plant/FlowerNumber), data = metDat)
plot(modAccMax)
exp(confint(modAccMax))


#########################################################################################
# Total pollen acceleration
pdf(file = "/Users/callinswitzer/Dropbox/ExperSummer2015/PolVel.pdf", width =8, height = 4)
for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #pollen
     pollen <- data.frame(x = dp$pt4_cam1_X, y= dp$pt4_cam1_Y)
     pollen <- pollen[ complete.cases(pollen), ]
     pollen$x <- pollen$x - min(pollen$x)
     pollen$y <- pollen$y - min(pollen$y)
     pollen <- pollen / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(pollen$x, spar = smspar) 
     ysmth <- smooth.spline(pollen$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     accMag <- with(acc, sqrt(xacc^2 + yacc^2))
     
     if(ii == 1){
          plot(y = velMag, x= ttme[1:(length(ttme)-1)], type = "l", 
               col = rgb(1,0,.3, 0.5),  xlim = c(-0.002, 0.1), ylim = c(-0, 5),
               xlab = "time (s)", ylab = "pollen velocity magnitude (m/s)")
     } 
     else lines(y = velMag, x= ttme[1:(length(ttme)-1)], col = rgb(1,0,0.3, 0.5))
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}
dev.off()


pdf("/Users/callinswitzer/Dropbox/ExperSummer2015/PolYAccel.pdf", width = 8, height = 4)
for(ii in 1:nrow(metDat)){
     ddfile <- paste(
          "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/", 
          metDat$digitizedFile[ii], sep = "")
     dp <- read.csv(ddfile)
     
     # calibration points
     pin <- data.frame(dp$pt1_cam1_X, dp$pt1_cam1_Y, dp$pt2_cam1_X, dp$pt2_cam1_Y)
     pin <- pin[complete.cases(pin), ]
     # to get to mm
     PixInPin <- (sqrt((pin$dp.pt1_cam1_X- pin$dp.pt2_cam1_X)^2 + 
                            (pin$dp.pt1_cam1_Y-pin$dp.pt2_cam1_Y)^2)) / metDat$CalSizeMM[ii] 
     
     #pollen
     pollen <- data.frame(x = dp$pt4_cam1_X, y= dp$pt4_cam1_Y)
     pollen <- pollen[ complete.cases(pollen), ]
     pollen$x <- pollen$x - min(pollen$x)
     pollen$y <- pollen$y - min(pollen$y)
     pollen <- pollen / PixInPin / 1000 # in meters
     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(pollen$x, spar = smspar) 
     ysmth <- smooth.spline(pollen$y, spar = smspar)
     
     # velocity
     vel <- data.frame(xvel = diff(xsmth$y), yvel = diff(ysmth$y))*5000
     acc <- with(vel, data.frame(xacc = diff(xvel), yacc = diff(yvel)))*5000
     
     velMag <- with(vel, sqrt(xvel^2 + yvel^2))
     
     indMax <- which(velMag == max(velMag))
     
     ttme <- -indMax:(length(velMag) - indMax) / 5000
     
     accMag <- with(acc, sqrt(xacc^2 + yacc^2))
     
     if(ii == 1){
          plot(y = acc$yacc, x= ttme[1:(length(ttme)-2)], type = "l", 
               col = rgb(0,0.3,1, 0.3),  xlim = c(-0.002, 0.1), ylim = c(-1500, 1000),
               xlab = "time (s)", ylab = "pollen y acceleration (m/s/s)", lwd = 1)
     } 
     else lines(y = acc$yacc, x= ttme[1:(length(ttme)-2)], col = rgb(0,0.3,1, 0.3), lwd = 1)
     #      plot(vel$yvel, type = "l")
     #      plot(acc$yacc, type = "l")
}

abline(h = -9.8, col = 'grey', lwd = 2)
legend('bottomright', legend = "Acc. Due to Gravity", lwd = 2, col = 'grey')
dev.off()
metDat
