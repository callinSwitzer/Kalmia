## Callin Switzer
## 29 Nov 2016
## Kalmia pollen and anther kinematics

# 1. Read in digitized files
# 2. Smooth digitized points, and impute
# 3. Use imputed points to calculate velocity and acceleration (normal and tangential)
# 4. 


## TODO:
# use cross-validation to decide smoothing parameters or plot noise/resoluation tradeoff
# or simply justify the choice -- by using visual inspection
# compute acceleration and velocity for different values of smoothing parameters
# idea -- show video with smoothed vs. unsmoothed points added -- background subtracted.



# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4", "signal", "reshape2", "viridis")

ipak(packages)


# read in metadata
dfile <- "/Users/callinswitzer/Dropbox/ExperSummer2015/LaurelsOnly.csv"
metDat <- read.csv(dfile)

metDat <- metDat[metDat$digitizedFile!= "", ]

# set constants:
fps <- 5000 # frames per second

ii = 11


# read in each .csv file for analysis
# make a list of data frames
digdirect <- "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/"

maxVals = data.frame()
for(kk in seq(from = 0.9, to = 0.05, by = -0.05)){

     newDF <- data.frame()
     
     for(ii in 1:nrow(metDat)){
          
          
          # ignore ii ==7, because the video started too late
          if(ii == 7) next
          
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
          
          # get frame where pollen starts and leaves
          antherPoll$polStart = 1:nrow(antherPoll) == metDat$framePollenStartedLeaving[ii]
          antherPoll$polEnd = 1:nrow(antherPoll) == metDat$framePollenReleaseComplete[ii]
          
          # gives only rows where either anth and pollen are complete
          antherPoll = antherPoll[ complete.cases(antherPoll[c('anthx')]) | 
                                        complete.cases(antherPoll[c('polx')]), ]
          
          # if x value starts to right of screen, reverse points,
          # so all x values start on the left part of the screen at 0
          if(lm(antherPoll[,1] ~ I(1:length( antherPoll[,1])))$coefficients[2] < 0 ){
               antherPoll$anthx <- metDat[ii,'vidWidth'] - antherPoll$anthx
               antherPoll$polx <- metDat[ii,'vidWidth'] - antherPoll$polx
          }
          
          
          # cbind data frame, to add smoothed columns
          antherPoll <- data.frame(cbind(antherPoll, antherPoll))
          
          # plot(x = antherPoll$anthx.1, y = antherPoll$anthy.1)
          # plot(antherPoll$anthx.1)
          
          # smooth with SG is based on the least-squares fitting of 
          # polynomials to segments of the data
          
          # other options include smoothing splines (tend to "cut the corners" of curves)
          # butterworth filters (inaccurate at endpoints)
          # Kernel smoothing
          
          x <- na.omit(antherPoll$anthy.1)
          xx <- c(x[round(length(x)/ 2):1], x, x[round(length(x)):round(length(x)/ 2)])
          want = c(rep(FALSE, round(length(x)/ 2)), rep(TRUE, length(x)), rep(FALSE, round(length(x)/2)))
          sg <- sgolayfilt(xx, p = 3, n = 11) # Savitzky-Golay filter
          # plot(xx[want], type="b", col = 'red', pch = 20)
          # points(sg[want], pch = 20, type = 'o') # smoothed SG data
          
          
          W = 0.99
          b1 <- butter(5, W, type = 'low')
          y1 <- filtfilt(b1, xx)
          
          # points(y1[want], pch=20, col='grey')
          
          
          
          
          
          # filter with Savitzky-Golay filter or Butterworth filter
          # degree = 3, frame size = 11 points
          foo = sapply(X = c("anthx.1", "anthy.1", "polx.1", "poly.1"), FUN = function(y){
               #sm1 <- sgolayfilt(na.omit(antherPoll[, x]), p = 3, n = 51) 
               
               # butterworth filter
               x <- na.omit(antherPoll[, y])
               xx <- c(x[round(length(x)/ 2):1], x, x[round(length(x)):round(length(x)/ 2)])
               want = c(rep(FALSE, round(length(x)/ 2)), rep(TRUE, length(x)), rep(FALSE, round(length(x)/2)))
               W = kk # sweet spot seems to be about 0.2
               b1 <- butter(5, W, type = 'low')
               y1 <- filtfilt(b1, xx)
               sm1 <- y1[want]
               antherPoll[, y][complete.cases(antherPoll[, y])] <<- sm1
          })
          
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
          
          
          
          # add columns to show velocity, based on smoothed, absolute position
          # velocity is in m/s
          bat = sapply(X = c("anthx.1.abs", "anthy.1.abs", "polx.1.abs", "poly.1.abs"), 
                       FUN = function(x){
                            newName = paste0(x, ".vel")
                            tmp <-  c(NaN, diff(antherPoll[,x])) * fps # add a NaN to beginning of data
                            antherPoll[,newName] <<- tmp 
                       })
          
          # calculate speed
          antherPoll$anthspeed = sqrt(antherPoll$anthx.1.abs.vel^2 + antherPoll$anthy.1.abs.vel^2)
          antherPoll$polspeed = sqrt(antherPoll$polx.1.abs.vel^2 + antherPoll$poly.1.abs.vel^2)
          
          # plot(antherPoll$anthspeed)
          
          ###########################################
          # pollen acceleration
          polVelocity = cbind(antherPoll$polx.1.abs.vel, antherPoll$poly.1.abs.vel)
          polSpeed = antherPoll$polspeed
          # plot(polSpeed)
          tme = antherPoll$tme
          polAccel = data.frame(rbind(c(NA, NA), apply(polVelocity, MARGIN = 2, FUN = diff))) * fps
          # par(mfrow =c(2,2))
          # plot(polAccel[,1], x = antherPoll$tme, type = 'l') # calculated
          # plot(polAccel[,2], x = antherPoll$tme, type = 'l')
          
          
          # unit tangent vector
          T_t = polVelocity / polSpeed
          
          DT = data.frame(rbind(c(NA, NA), apply(T_t, MARGIN = 2, FUN = diff))) * fps
          NormDT = sqrt(DT[,1]^2 + DT[,2]^2)
          Curvature = NormDT / polSpeed
          
          # compute a_N (normal acceleration) and a_T (tangential acceleration)
          # a_T = ds/dt
          a_T =  c(NA, diff(polSpeed) * fps)
          N_t = data.frame(t(sapply(1:nrow(DT), FUN = function(x) unlist(DT[x, ] / NormDT[x]))))
          # plot(a_T, type = "l", ylim = c(-3000, 3000))
          
          # a_N = speed^2 * curvature
          a_N = polSpeed^2 * Curvature
          
          
          # check total accel by adding normal and tangential accelerations
          # a_total = a_T * T_t + a_N * N_t
          a_total = as.data.frame(t(sapply(X = 1:nrow(polAccel), FUN  = function(x) a_T[x] * T_t[x, ] + a_N[x] * N_t[x,] )))
          # plot(a_total) # includes both x and y coordinates
          # plot(polAccel)
          
          # par(mfrow = c(2,2))
          # plot(unlist(a_total[,1]))
          # plot(unlist(a_total[,2]))
          # plot(polAccel[,1])
          # plot(polAccel[,2])
          
          
          # plot(a_N)
          # plot(a_T)
          # 
          a_T_Pol = a_T
          # plot(a_T, x = tme)
          # plot(a_N, x = tme)
          
          
          # calculate magnitude of acceleration, using two methods
          # 1. Normal and tangential acceleration
          a_mag1 = sqrt(a_T^2 + a_N^2)
          # plot(a_mag1)
          
          amag2 = sqrt(polAccel[,1]^2 + polAccel[,2]^2)
          # plot(amag2, type = 'l')
          # plot(polVelocity[,1])
          
          
          ########################################
          
          ###########################################
          # anther acceleration
          anthVelocity = cbind(antherPoll$anthx.1.abs.vel, antherPoll$anthy.1.abs.vel)
          anthSpeed = antherPoll$anthspeed
          # plot(anthSpeed)
          tme = antherPoll$tme
          anthAccel = data.frame(rbind(c(NA, NA), apply(anthVelocity, MARGIN = 2, FUN = diff))) * fps
          # par(mfrow =c(2,2))
          # # plot(anthAccel[,1], x = antherPoll$tme, type = 'l') # calculated
          # plot(anthAccel[,2], x = antherPoll$tme, type = 'l')
          
          
          # unit tangent vector
          T_t = anthVelocity / anthSpeed
          
          DT = data.frame(rbind(c(NA, NA), apply(T_t, MARGIN = 2, FUN = diff))) * fps
          NormDT = sqrt(DT[,1]^2 + DT[,2]^2)
          Curvature = NormDT / anthSpeed
          
          # compute a_N (normal acceleration) and a_T (tangential acceleration)
          # a_T = ds/dt
          a_T =  c(NA, diff(anthSpeed) * fps)
          a_T_anth = a_T
          N_t = data.frame(t(sapply(1:nrow(DT), FUN = function(x) unlist(DT[x, ] / NormDT[x]))))
          # plot(a_T, type = "l", ylim = c(-3000, 3000))
          
          # a_N = speed^2 * curvature
          a_N = anthSpeed^2 * Curvature
          
          
          # check total accel by adding normal and tangential accelerations
          # a_total = a_T * T_t + a_N * N_t
          a_total = as.data.frame(t(sapply(X = 1:nrow(anthAccel), FUN  = function(x) a_T[x] * T_t[x, ] + a_N[x] * N_t[x,] )))
          # plot(a_total) # includes both x and y coordinates
          # plot(anthAccel)
          
          # par(mfrow = c(2,2))
          # plot(unlist(a_total[,1]))
          # plot(unlist(a_total[,2]))
          # plot(anthAccel[,1])
          # plot(anthAccel[,2])
          
          
          # plot(a_N)
          # plot(a_T)
          
          # par(mfrow = c(2,1))
          # plot(a_T, x = tme, type = 'l')
          # max(a_T, na.rm = TRUE)
          # which.max(a_T)
          
          # plot(a_T_Pol, x = tme, type = 'l')
          # max(a_T_Pol, na.rm = TRUE)
          # which.max(a_T_Pol)
          # 
          tmeRoll <- seq(from = -which(antherPoll$polStart) + 1, length.out = length(tme)) / fps
          
          dfi <- data.frame(anthSpeed, polSpeed, a_T_anth, a_T_Pol, tme, 
                            trial = metDat$VideoName[ii], 
                            tmeStart = antherPoll$polStart,
                            tmeEnd = antherPoll$polEnd, 
                            centeredTime = tmeRoll)
          
          newDF <- rbind(newDF,dfi) 
          print(ii)
          
     }


antherPoll$frame <- 1:nrow(antherPoll)
ggplot(na.omit(antherPoll)) + 
     geom_point(aes(x = anthx, y = anthy), colour = "grey", alpha = 0.3) + 
     geom_path(aes(x = anthx, y = anthy), color = "grey", alpha = 0.3) + 
     geom_point(aes(x = anthx.1, y = anthy.1, colour = frame)) + 
     geom_path(aes(x = anthx.1, y = anthy.1, color = frame)) + 
     scale_color_viridis()

ggplot((antherPoll)) + 
     geom_point(aes(x = polx, y = poly), colour = "grey", alpha = 0.3) + 
     geom_path(aes(x = polx, y = poly), color = "grey", alpha = 0.3) + 
     geom_point(aes(x = polx.1, y = poly.1, colour = frame)) + 
     geom_path(aes(x = polx.1, y = poly.1, color = frame)) + 
     scale_color_viridis() + 
     coord_fixed(ratio = 1)

plot(antherPoll$anthx, antherPoll$anthy)
points(antherPoll$anthx.1, antherPoll$anthy.1, type = 'b', pch = 20)
     



theme_set(theme_classic())

savePath = "/Users/callinswitzer/Dropbox/ExperSummer2015/Kalmia2015FiguresAndData/"

# anther speed
ggplot(newDF, aes(x = centeredTime, y = anthSpeed, group = trial)) + 
     geom_line(alpha = 0.5) +
     xlim(c(-0.01, 0.02)) + 
     ylim(c(0,6)) + 
     labs(x = "Time (s)", y = "Anther speed (m/s)")
ggsave(paste0(savePath, "antherSpeed", "filt_", kk, ".pdf"), width = 5, height = 4)


# pollen speed
ggplot(newDF, aes(x = centeredTime, y = polSpeed, group = trial)) + 
     geom_line(alpha = 0.5) + 
     xlim(c(-0.01, 0.02)) +
     ylim(c(0,6)) +  
     labs(x = "Time (s)", y = "Pollen speed (m/s)")
ggsave(paste0(savePath, "pollenSpeed", "filt_", kk, ".pdf"), width = 5, height = 4)

# anther tangential acceleration
ggplot(newDF, aes(x = centeredTime, y = a_T_anth, group = trial)) + 
     geom_line(alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Anther tangential acceleration  (m/s/s)")
ggsave(paste0(savePath, "antherTangAccel", "filt_", kk, ".pdf"), width = 5, height = 4)

# pollen tangential acceleration
# anther tangential acceleration
ggplot(newDF, aes(x = centeredTime, y = a_T_Pol, group = trial)) + 
     geom_line(alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Pollen tangential acceleration  (m/s/s)")
ggsave(paste0(savePath, "PollenTangAccel", "filt_", kk, ".pdf"), width = 5, height = 4)

# find max for each measurement for each trial

# anther speed
mmx <- as.data.frame(t(sapply(unique(as.character(newDF$trial)), FUN = function(x){
     tmp <- newDF[newDF$trial == x, ]
     return (unlist(tmp[which.max(tmp$anthSpeed),]))
})))
mmx$trial <- row.names(mmx)


ggplot() + 
     geom_line(data = newDF, aes(x = centeredTime, y = anthSpeed, group = as.factor(trial)), alpha = 0.5) +
     xlim(c(-0.01, 0.02)) + 
     ylim(c(0,6)) + 
     labs(x = "Time (s)", y = "Anther speed (m/s)") + 
     geom_point(data = mmx, aes(x = centeredTime, y = anthSpeed), color = 'red', alpha = 0.5) + 
     theme(legend.position = "none") 
#+  facet_wrap(~ trial)
ggsave(paste0(savePath, "antherSpeedMax", "filt_", kk, ".pdf"), width = 5, height = 4)


# pollen speed
mmp <- as.data.frame(t(sapply(unique(as.character(newDF$trial)), FUN = function(x){
     tmp <- newDF[newDF$trial == x, ]
     tmp <- tmp[abs(tmp$centeredTime) < 0.01, ]
     return (unlist(tmp[which.max(tmp$polSpeed),]))
})))

mmp$trial <- row.names(mmp)

# pollen speed
ggplot() + 
     geom_line(data = newDF, aes(x = centeredTime, y = polSpeed, group = trial), alpha = 0.5) + 
     xlim(c(-0.01, 0.02)) + 
     ylim(c(0,6)) + 
     labs(x = "Time (s)", y = "Pollen speed (m/s)") + 
     geom_point(data = mmp, aes(x = centeredTime, y = polSpeed), color = 'red', alpha = 0.5) + 
     theme(legend.position = "none") 
ggsave(paste0(savePath, "pollenSpeedMax", "filt_", kk, ".pdf"), width = 5, height = 4)
     


# anther acceleration
mma <- as.data.frame(t(sapply(unique(as.character(newDF$trial)), FUN = function(x){
     tmp <- newDF[newDF$trial == x, ]
     
     # get only points that are within 0.05 seconds of the centered time
     # to ignore the anthers hitting the other side of the flower
     tmp <- tmp[abs(tmp$centeredTime) < 0.005, ]
     return (unlist(tmp[which.max(tmp$a_T_anth),]))
})))

mma$trial <- row.names(mma)

ggplot() + 
     geom_line(data = newDF, aes(x = centeredTime, y = a_T_anth, group = trial), alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Anther tangential acceleration  (m/s/s)") + 
     geom_point(data = mma, aes(x = centeredTime, y = a_T_anth), color = 'red', alpha = 0.5)
ggsave(paste0(savePath, "antherTangAccelMax", "filt_", kk, ".pdf"), width = 5, height = 4)


# pollen acceleration
mmpp <- as.data.frame(t(sapply(unique(as.character(newDF$trial)), FUN = function(x){
     tmp <- newDF[newDF$trial == x, ]
     tmp <- tmp[abs(tmp$centeredTime) < 0.007, ]
     return (unlist(tmp[which.max(tmp$a_T_Pol),]))
})))

mmpp$trial <- row.names(mmpp)

ggplot() + 
     geom_line(data = newDF, aes(x = centeredTime, y = a_T_Pol, group = trial), alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Pollen tangential acceleration  (m/s/s)") + 
     geom_point(data = mmpp, aes(x = centeredTime, y = a_T_Pol), color = 'red', alpha = 0.5)

ggsave(paste0(savePath, "pollenTangAccelMax", "filt_", kk, ".pdf"), width = 5, height = 4)


# estimate ranges for acceleration, and speed
md  = merge(x = mmx[, c('trial', 'anthSpeed')], metDat, by.x = "trial", by.y = 
                 "VideoName")
md = merge(x = mmp[, c('trial', 'polSpeed')], md, by = "trial")

md = merge(x = mmpp[, c('trial', 'a_T_Pol')], md, by = "trial")
md = merge(x = mma[, c('trial', 'a_T_anth')], md, by = "trial")


#LMER

modVelMaxAnth <- lmer(formula = anthSpeed ~ (1|plant/FlowerNumber), data = md)
summary(modVelMaxAnth)
#confint(modVelMaxAnth)


modVelMaxPol <- lmer(formula = polSpeed ~ (1|plant/FlowerNumber), data = md)
summary(modVelMaxPol)
#confint(modVelMaxPol)


modAccMaxPol <- lmer(formula = a_T_Pol ~ (1|plant/FlowerNumber), data = md)
summary(modAccMaxPol)
#confint(modAccMaxPol)


modAccMaxAnth <- lmer(formula = a_T_anth ~ (1|plant/FlowerNumber), data = md)
summary(modAccMaxAnth)
#confint(modAccMaxAnth)



newRow = c(summary(modVelMaxAnth)$coef[1], summary(modVelMaxPol)$coef[1], 
           summary(modAccMaxPol)$coef[1], summary(modAccMaxAnth)$coef[1])

maxVals <- rbind(maxVals, newRow)
}


colnames(maxVals) <- c("Max Veloc Anth (m/s)", "Max Veloc Pol (m/s)", "Max Acc Pol (m/s)", "Max Acc Anth (m/s)")
maxVals$filtParameter <- seq(from = 0.9, to = 0.05, by = -0.05)
maxVals[, is.na(colnames(maxVals)) ] <- NULL


maxV_long <- melt(maxVals, id.vars = "filtParameter")
maxV_long$acc = as.numeric(sapply(maxV_long$variable, function(x) grep(pattern = "Acc", as.character(x)) == 1))
library(plyr)
maxV_long$acc <- mapvalues(as.character(maxV_long$acc), from = c("1", NA), to = c("Acceleration", "Velocity"))


ggplot(maxV_long, aes(x = filtParameter, y = value, color = variable)) + 
     geom_line(size = 2) + 
     facet_wrap(~acc, scales = 'free') + 
     labs(x = "filter parameter", y = "Value") + 
     scale_color_viridis(discrete = TRUE, name = "Measured Variable")

