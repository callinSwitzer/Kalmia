## Callin Switzer
## 29 Nov 2016
## 6 Dec 2016 Update -- re-digitized pollen and anthers very carefully
## may not need any smoothing
## 8 Dec Update: used cross validation for decide smoothing parameter for smoothing spline
## 8 Feb 2017 Update: Cleaned up code for submission as supplemental info



## Kalmia pollen and anther kinematics
# 1. Read in digitized files
# 2. Smooth digitized points, using smoothing spline with tuning parameter from cross validation
# 3. Use smoothed points to calculate velocity and acceleration (normal and tangential)



# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4", "signal")

ipak(packages)


# read in metadata
dfile <- "/Users/callinswitzer/Dropbox/ExperSummer2015/LaurelsOnly.csv"
metDat <- read.csv(dfile)

metDat <- metDat[metDat$redigitizedFile!= "", ]

# set constants:
fps <- 5000 # frames per second

ii = 11


# read in each .csv file for analysis
# make a list of data frames
digdirect <- "/Users/callinswitzer/Dropbox/ExperSummer2015/CleanKalmiaDigitized/"


newDF <- data.frame()

for(ii in 1:nrow(metDat)){
     
     ddfile <- paste0(digdirect, metDat$redigitizedFile[ii])
     #ddfile <- paste0(digdirect, "20150611_120046xypts.csv")
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
     
     # smooth with a smoothing spline, spar = 0.29
     foo = sapply(X = c("anthx.1", "anthy.1", "polx.1", "poly.1"), FUN = function(y){
          yy = na.omit(antherPoll[, y])
          xx = 1:length(yy)
          
          sm1 <- smooth.spline(x = xx, y = yy, spar = 0.29)
          antherPoll[, y][complete.cases(antherPoll[, y])] <<- sm1$y
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
     
     # plot(antherPoll$anthspeed, type = 'l')
     
     # plot(antherPoll$anthx.1, antherPoll$anthy.1)
     # points(antherPoll$anthx, antherPoll$anthy, pch = 20)
     # plot(antherPoll$anthx)
     # plot(antherPoll$anthy)
     
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

theme_set(theme_classic())

savePath = "/Users/callinswitzer/Dropbox/ExperSummer2015/Kalmia2015FiguresAndData/"

# anther speed
ggplot(newDF, aes(x = centeredTime, y = anthSpeed, group = trial)) + 
     geom_line(alpha = 0.5) +
     xlim(c(-0.01, 0.02)) + 
     ylim(c(0,6)) + 
     labs(x = "Time (s)", y = "Anther speed (m/s)")
ggsave(paste0(savePath, "antherSpeed01_CVSmoothSpline.pdf"), width = 5, height = 4)


# pollen speed
ggplot(newDF, aes(x = centeredTime, y = polSpeed, group = trial)) + 
     geom_line(alpha = 0.5) + 
     xlim(c(-0.01, 0.02)) +
     ylim(c(0,6)) +  
     labs(x = "Time (s)", y = "Pollen speed (m/s)")
ggsave(paste0(savePath, "pollenSpeed01_CVSmoothSpline.pdf"), width = 5, height = 4)

# anther tangential acceleration
ggplot(newDF, aes(x = centeredTime, y = a_T_anth, group = trial)) + 
     geom_line(alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Anther tangential acceleration  (m/s/s)")
ggsave(paste0(savePath, "antherTangAccel01_CVSmoothSpline.pdf"), width = 5, height = 4)

# pollen tangential acceleration
# anther tangential acceleration
ggplot(newDF, aes(x = centeredTime, y = a_T_Pol, group = trial)) + 
     geom_line(alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Pollen tangential acceleration  (m/s/s)")
ggsave(paste0(savePath, "PollenTangAccel01_CVSmoothSpline.pdf"), width = 5, height = 4)

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
ggsave(paste0(savePath, "antherSpeedMax01_CVSmoothSpline.pdf"), width = 5, height = 4)


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
ggsave(paste0(savePath, "pollenSpeedMax01_CVSmoothSpline.pdf"), width = 5, height = 4)



# anther acceleration
mma <- as.data.frame(t(sapply(unique(as.character(newDF$trial)), FUN = function(x){
     tmp <- newDF[newDF$trial == x, ]
     
     # get only points that are within 0.05 seconds of the centered time
     # to ignore the anthers hitting the other side of the flower
     #tmp <- tmp[abs(tmp$centeredTime) < 0.002, ]
     return (unlist(tmp[which.max(tmp$a_T_anth),]))
})))

mma$trial <- row.names(mma)

ggplot() + 
     geom_line(data = newDF, aes(x = centeredTime, y = a_T_anth, group = trial), alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Anther tangential acceleration  (m/s/s)") + 
     geom_point(data = mma, aes(x = centeredTime, y = a_T_anth), color = 'red', alpha = 0.5)

ggsave(paste0(savePath, "antherTangAccelMax01_CVSmoothSpline.pdf"), width = 5, height = 4)


# pollen acceleration
mmpp <- as.data.frame(t(sapply(unique(as.character(newDF$trial)), FUN = function(x){
     tmp <- newDF[newDF$trial == x, ]
     #tmp <- tmp[abs(tmp$centeredTime) < 0.003, ]
     return (unlist(tmp[which.max(tmp$a_T_Pol),]))
})))

mmpp$trial <- row.names(mmpp)

ggplot() + 
     geom_line(data = newDF, aes(x = centeredTime, y = a_T_Pol, group = trial), alpha = 0.5) + 
     #ylim(c(-2500, 4000)) +
     xlim(c(-0.01, 0.02)) + 
     labs(x = "Time (s)", y = "Pollen tangential acceleration  (m/s/s)") + 
     geom_point(data = mmpp, aes(x = centeredTime, y = a_T_Pol), color = 'red', alpha = 0.5)

ggsave(paste0(savePath, "pollenTangAccelMax01_CVSmoothSpline.pdf"), width = 5, height = 4)


# estimate ranges for acceleration, and speed
md  = merge(x = mmx[, c('trial', 'anthSpeed')], metDat, by.x = "trial", by.y = 
                 "VideoName")
md = merge(x = mmp[, c('trial', 'polSpeed')], md, by = "trial")

md = merge(x = mmpp[, c('trial', 'a_T_Pol')], md, by = "trial")
md = merge(x = mma[, c('trial', 'a_T_anth')], md, by = "trial")


#LMER
hist(md$anthSpeed) #dist is fine
modVelMaxAnth <- lmer(formula = anthSpeed ~  (1|plant/FlowerNumber), data = md)
summary(modVelMaxAnth)
confint(modVelMaxAnth)

# diagnostics
plot(modVelMaxAnth)

# QQPlot for group-level effects
qqnorm(ranef(modVelMaxAnth)$plant[[1]], main="Normal Q-Q plot for random effects")
qqline(ranef(modVelMaxAnth)$plant[[1]])

# QQPlot for group-level effects
qqnorm(ranef(modVelMaxAnth)$FlowerNumber[[1]], main="Normal Q-Q plot for random effects")
qqline(ranef(modVelMaxAnth)$FlowerNumber[[1]]) # 



modVelMaxPol <- lmer(formula = polSpeed ~ (1|plant/FlowerNumber), data = md)
summary(modVelMaxPol)
confint(modVelMaxPol)
plot(modVelMaxPol)


# log transformed, b/c distribution is skewed.
hist(md$a_T_Pol)
hist(log(md$a_T_Pol)) #better
modAccMaxPol <- lmer(formula = log(a_T_Pol) ~  (1|plant/FlowerNumber), data = md)
summary(modAccMaxPol)
exp(confint(modAccMaxPol))
plot(modAccMaxPol)

modAccMaxAnth <- lmer(formula = log(a_T_anth) ~ (1|plant/FlowerNumber), data = md)
summary(modAccMaxAnth)
exp(confint(modAccMaxAnth))
plot(modAccMaxAnth)
# diagnostics
plot(modAccMaxAnth)

# QQPlot for group-level effects
qqnorm(ranef(modAccMaxAnth)$plant[[1]], main="Normal Q-Q plot for random effects")
qqline(ranef(modAccMaxAnth)$plant[[1]])

# QQPlot for group-level effects
qqnorm(ranef(modAccMaxAnth)$FlowerNumber[[1]], main="Normal Q-Q plot for random effects")
qqline(ranef(modAccMaxAnth)$FlowerNumber[[1]]) # 

