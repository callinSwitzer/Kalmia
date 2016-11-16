## Laurel Digitizations
# Callin Switzer
# June 21, 2015
# Updated 16 Nov 2016 - redone to make code better
# Calculate velocity, position, and acceleration, of anther and pollen release 
# from mountain laurels from Arboretum.


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

#for(ii in 1:length(metDat$digitizedFile)){
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
     
     # get anther locations
     anther <- data.frame(x = dp$pt3_cam1_X, y= dp$pt3_cam1_Y)
     anther = anther[ complete.cases(anther), ]
     


#




















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
     par(pty="s")
     if(ii == 1) plot(anther$x, anther$y, ylab = "y position (meters)", xlab = "X position (meters)", 
                      xlim = c(0, 0.01), ylim = c(0, 0.01), pch = ".")
     points(anther$x, anther$y, pch = ".")
     tme <- 1:length(anther$x) / fps # time
     

     
     # Smooth
     smspar <- 0.5
     xsmth <- smooth.spline(anther$x, spar = smspar) 
     ysmth <- smooth.spline(anther$y, spar = smspar)
     #?smooth.spline
     
     lines(xsmth$y, ysmth$y, type = "l", lwd = 3, col = rgb(0, 0, 0, 0.5)) 
     
     Sys.sleep(0.2)
     
}


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
