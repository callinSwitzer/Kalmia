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



# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4", "signal", "mgcv")

ipak(packages)


# read in metadata
dfile <- "/Users/callinswitzer/Dropbox/ExperSummer2015/LaurelsOnly.csv"
metDat <- read.csv(dfile)

metDat <- metDat[metDat$digitizedFile!= "", ]

# set constants:
fps <- 5000 # frames per second

ii = 16

# read in each .csv file for analysis
# make a list of data frames
digdirect <- "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/"
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


# cbind data frame, to add smoothed columns
antherPoll <- data.frame(cbind(antherPoll, antherPoll))

plot(x = antherPoll$anthx.1, y = antherPoll$anthy.1)
plot(antherPoll$anthx.1)







# smooth with SG is based on the least-squares fitting of 
# polynomials to segments of the data

# other options include smoothing splines (tend to "cut the corners" of curves)
# butterworth filters (inaccurate at endpoints)
# Kernel smoothing

x <- na.omit(antherPoll$poly.1)
sg <- sgolayfilt(x, p = 3, n = 11) # Savitzky-Golay filter
plot(x, type="b", col = 'red', pch = 20)
lines(sg, pch = 20, type = 'l') # smoothed SG data
# interpolation
inn = interp1(x = 1:length(x), xi = seq(1, length(x), length.out = 500), y = sg, method = "cubic")
lines(inn, col = 'pink', x = seq(1, length(x), length.out = 500), pch = 20, type = 'b')

# look at residuals
plot(x - sg, type = 'b')
abline(h = 0, lty = 2)

# calculate run lengths (bias)
rrll = rle(x - sg > 0)
mean(rrll$lengths)
hist(log(rrll$lengths))


# evaluate bias
# want to minimize long runs with high distance from line 
## HERE -- add up distance from line in runs
aa <- numeric()
bb = numeric()
cc = numeric()
dd= numeric()
xx = seq(5, 91, by = 2)
for(ii in xx){
     x <- na.omit(antherPoll$poly.1)
     sg <- sgolayfilt(x, p = 3, n = ii) # Savitzky-Golay filter
     # calculate run lengths (bias)
     rrll = rle(x - sg > 0)
     aa <- c(aa, mean(log(rrll$lengths)))
     
     bb <- c(bb, sd(sg - smooth.spline(x)$y))
     cc <- c(cc, sqrt(1 / length(x) * sum((x - sg)^2))) # Root mean squared error
     dd <- c(dd, mad(x - sg))
     
}
plot(y = aa, xx)
abline(v = 11)
plot(y = bb, xx)
abline(v = 11)
plot(y = cc, xx)
abline(v = 11)
plot(y = dd, xx)
abline(v = 11)


cumsum(rrll$lengths)


# interpolation of raw data
inn = interp1(x = 1:length(x), xi = seq(1, length(x), length.out = 500), y = x, method = "cubic")
lines(inn, col = 'grey', x = seq(1, length(x), length.out = 500), pch = 20, type = 'b')


# filter with Savitzky-Golay filter
# degree = 3, frame size = 11 points
foo = sapply(X = c("anthx.1", "anthy.1", "polx.1", "poly.1"), FUN = function(x){
     sm1 <- sgolayfilt(na.omit(antherPoll[, x]), p = 3, n = 11) 
     antherPoll[, x][complete.cases(antherPoll[, x])] <<- sm1
})
plot(x = antherPoll$anthx.1, y = antherPoll$anthy.1)
plot(x = antherPoll$anthx, y = antherPoll$anthy)

plot( antherPoll$anthx)
points(antherPoll$anthx.1, pch = 20)
plot( antherPoll$anthy)
points(antherPoll$anthy.1, pch = 20)


# recently smoothed data
plot(antherPoll$poly.1, col = 'blue')



# interpolation
inn = interp1(x = 1:length(x), xi = seq(1, length(x), by = 0.1), y = sg, method = "cubic")
lines(inn, col = 'pink', x = seq(1, length(x), by = 0.1), pch = 20, type = 'b')
length(x)
length(inn)


# look at residuals
plot(x - sg, type = 'b')
plot(x - na.omit(antherPoll$poly.1), type = 'b')
abline(h = 0, lty = 2)

plot(y = x, x = 1:length(x))
