# acceleration calculating examples
# radial coordinates trying to calculate total acceleration


# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4", "plotrix")

ipak(packages); rm(packages)


# fit circle to data
# function from here: 
# http://stackoverflow.com/questions/27169122/how-to-find-a-best-fit-circle-ellipse-using-r
fitSS <- function(xy,
                  a0=mean(xy[,1]),
                  b0=mean(xy[,2]),
                  r0 = mean(sqrt((xy[,1]-a0)^2 + (xy[,2]-b0)^2)),
                  ...){
     SS <- function(abr){
          sum((abr[3] - sqrt((xy[,1]-abr[1])^2 + (xy[,2]-abr[2])^2))^2)
     }
     optim(c(a0,b0,r0), SS, ...)
}

fit.ellipse <- function (x, y = NULL) {
     # from:
     # http://r.789695.n4.nabble.com/Fitting-a-half-ellipse-curve-tp2719037p2720560.html
     #
     # Least squares fitting of an ellipse to point data
     # using the algorithm described in: 
     #   Radim Halir & Jan Flusser. 1998. 
     #   Numerically stable direct least squares fitting of ellipses. 
     #   Proceedings of the 6th International Conference in Central Europe 
     #   on Computer Graphics and Visualization. WSCG '98, p. 125-132 
     #
     # Adapted from the original Matlab code by Michael Bedward (2010)
     # michael.bedward@gmail.com
     #
     # Subsequently improved by John Minter (2012)
     # 
     # Arguments: 
     # x, y - x and y coordinates of the data points.
     #        If a single arg is provided it is assumed to be a
     #        two column matrix.
     #
     # Returns a list with the following elements: 
     #
     # coef - coefficients of the ellipse as described by the general 
     #        quadratic:  ax^2 + bxy + cy^2 + dx + ey + f = 0 
     #
     # center - center x and y
     #
     # major - major semi-axis length
     #
     # minor - minor semi-axis length
     #
     EPS <- 1.0e-8 
     dat <- xy.coords(x, y) 
     
     D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y) 
     D2 <- cbind(dat$x, dat$y, 1) 
     S1 <- t(D1) %*% D1 
     S2 <- t(D1) %*% D2 
     S3 <- t(D2) %*% D2 
     T <- -solve(S3) %*% t(S2) 
     M <- S1 + S2 %*% T 
     M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2) 
     evec <- eigen(M)$vec 
     cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2 
     a1 <- evec[, which(cond > 0)] 
     f <- c(a1, T %*% a1) 
     names(f) <- letters[1:6] 
     
     # calculate the center and lengths of the semi-axes 
     #
     # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2288654/
     # J. R. Minter
     # for the center, linear algebra to the rescue
     # center is the solution to the pair of equations
     # 2ax +  by + d = 0
     # bx  + 2cy + e = 0
     # or
     # | 2a   b |   |x|   |-d|
     # |  b  2c | * |y| = |-e|
     # or
     # A x = b
     # or
     # x = Ainv b
     # or
     # x = solve(A) %*% b
     A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
     b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
     soln <- solve(A) %*% b
     
     b2 <- f[2]^2 / 4
     
     center <- c(soln[1], soln[2]) 
     names(center) <- c("x", "y") 
     
     num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6]) 
     den1 <- (b2 - f[1]*f[3]) 
     den2 <- sqrt((f[1] - f[3])^2 + 4*b2) 
     den3 <- f[1] + f[3] 
     
     semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) )) 
     
     # calculate the angle of rotation 
     term <- (f[1] - f[3]) / f[2] 
     angle <- atan(1 / term) / 2 
     #angle <- atan2(f[2],f[1]-f[3]) / 2 # fixed angle calculation (sometimes gives wrong sign)
     list(coef=f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle)) 
}

get.ellipse <- function( fit, n=360 ) 
{
     # Calculate points on an ellipse described by 
     # the fit argument as returned by fit.ellipse 
     # 
     # n is the number of points to render 
     
     tt <- seq(0, 2*pi, length=n) 
     sa <- sin(fit$angle) 
     ca <- cos(fit$angle) 
     ct <- cos(tt) 
     st <- sin(tt) 
     
     x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa 
     y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca 
     
     cbind(x=x, y=y) 
}

efit <- fit.ellipse()

# read in metadata
dfile <- "/Users/callinswitzer/Dropbox/ExperSummer2015/LaurelsOnly.csv"
metDat <- read.csv(dfile)

metDat <- metDat[metDat$digitizedFile!= "", ]

# set constants:
fps <- 5000 # frames per second

# read in each .csv file for analysis
# make a list of data frames
digdirect <- "/Users/callinswitzer/Dropbox/ExperSummer2015/AllLaurelsDigitizations/"


########################## EXAMPLE 1 #############################
## example points around a circle
## constant tangential velocity
## tangential acceleration is 0, but normal acceleration is constant
pi = 3.14159
fps = 1000
circlePoints = function(r, n ){
     theta  = (pi*2 * (0:n)/n)
     return(data.frame(x = r * cos(theta), y = r * sin(theta)))
}

foo <- circlePoints(1, 1000)
newD <- foo[foo[,2] > 0, ]
newD <- newD[nrow(newD):1, ]

# different example data
# x position is changing in even intervals
# y is not -- non constant acceleration
# newD <- data.frame(x = 1:3000 / 1000, y = sin(1:3000/1000))
newD <- within(newD, {
     x = x - min(x)
     y = y - min(y)
})


plot(newD, asp = 1)
plot(newD$x)
plot(newD$y)
tme <- (1:length(newD$x)) / fps


# tangential displacement
# total distance traveled (by integrating the distance traveled in x and y, combined)
dd1 <- sqrt(diff(newD$x)^2 + 
                 diff(newD$y)^2)
dd1[is.na(dd1)] <- 0
dist <- c(0, cumsum(dd1))
# distance traveled
plot(y = dist, x = tme)
#points(newD$y, col = 'red')
#points(newD$x, col = 'blue')
points(c(0, cumsum(diff(newD$x))), x = tme, col = 'blue')
points(c(0, cumsum(abs(diff(newD$y)))), x = tme, col = 'red')
legend("topleft", legend = c("total distance", "y distance", "x distance"), col = c('black', 'red', 'blue'), pch = 1)
tail(dist)


# total displacement (from origin)
dd1 <- sqrt((newD$x - newD$x[1])^2 + 
                 (newD$y - newD$y[1])^2)
dd1[is.na(dd1)] <- 0
# distance traveled
plot(dd1, x = tme)
points(newD$y, x = tme,  col = 'red')
points(newD$x, x = tme, col = 'blue')
legend("topleft", legend = c("total displacement", "y disp.", "x disp"), col = c('black', 'red', 'blue'), pch = 1)


# tangential or total? speed (goes with distance, not displacement)
#calculate one way
tv <- c(NA, diff(dist)*fps)
plot(round(tv, 5))

# calculate another way (same values)
tx <- c(NA, diff(newD$x)*fps)
ty <- c(NA, diff(newD$y)*fps)
tv1 <- sqrt(tx^2 + ty^2)
plot(round(tv1, 5))

# all equal
all(round(tv1, 5) == round(tv, 5), na.rm = TRUE)


x = tme
par(mfrow = c(2,1), mai = c(0,1, 0,0))
plot(x, dist)
curve(1.45*x, add = TRUE, col = 'green')
plot(round(tv, 5), col = 'red', ylab = "velocity") # should be constant



# tangential acceleration -- should always be 0, because
# speed is constant
ta <- c(NA, diff(tv)*fps)
plot(round(ta, 5))
tail(ta)

# calculate total accel, and see if it's the same
ax <- c(NA, diff(tx)*fps)
ay <- c(NA, diff(ty)*fps)
plot(ax)
plot(ay)
ta1 <- sqrt(ax^2 + ay^2)
plot(ta1) # total accel


# normal acceleration -- should be constant, and not zero, for circle
# a = v^2 / r
# calculate normal acceleration with changing radius
circ3Pt <- function(xy){
     # returns center and radius of circle, given three points
     if(nrow(xy) != 3){
          return(c(NA, NA, NA))
     }
     x = xy[, 1]
     y = xy[, 2]
     
     mr <- (y[2] - y[1]) / (x[2] - x[1])
     mt <- (y[3] - y[2]) / (x[3] - x[2])
     
     xc <- (mr * mt * (y[3] - y[1]) + mr * (x[2] + x[3]) - mt* (x[1] + x[2])) / (2 * (mr - mt))
     yc <- (-1 / mr) * (xc - (x[1] + x[2]) / 2) + (y[1] + y[2]) / 2
     
     rad = sqrt((x[2] - xc)^2 + (y[2] - yc)^2)
     draw.circle(xc, yc,rad)
     return(c(xc, yc, rad))
     
     
}
dev.off()
plot(newD$x, newD$y, ylim = c(-5, 1), asp = 1)


# now right
centRad = t(sapply(X = 1:length(newD$x), FUN = function(x){
     tryCatch(circ3Pt(newD[(x-2):(x), ]), 
              error = function(e){return(c(NA, NA, NA))})
}))

# calculate total acceleration with changing radius
# normal acceleration -- should be constant, and not zero
# a = v^2 / r
aN <- ta^2 / centRad[,3]
plot(round(aN, 5))
plot(ta1)

plot(ta)
points(aN)
plot(aN)
tail(aN)
tta <- sqrt(ta1^2 + aN^2)
tail(tta)
tail(ta1)

plot(tta)
points(ta1)



plot(tv)
plot(round(ta1, 5))
points(round(tta, 5))
tail(ta1) 
tail(tta) # same, for circle










# example data
# x position is changing in even intervals
# y is not -- non constant acceleration
df1 <- data.frame(x = 1:30 / 10, y = sin(1:30/10))
# position
df1 <- within(df1, {
     x = x - min(x)
     y = y - min(y)
})

plot(df1)

# total distance
dd1 <- sqrt(diff(df1$x)^2 + 
                 diff(df1$y)^2)
dd1[is.na(dd1)] <- 0
dist <- c(0, cumsum(dd1))
# distance traveled
plot(y = dist, x = df1$x)
points(df1)



plot(df1$x)
plot(df1$y)




ii = 0

{
     ii = ii + 1
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
     
     #play video
     plot(c(antherPoll$polx,antherPoll$anthx),
          c(antherPoll$poly, antherPoll$anthy),
          type = 'n', asp = 1)
     
     # slow method, but can make a movie, if needed
     for(jj in 1:nrow(antherPoll)){
          points(antherPoll$polx[jj], antherPoll$poly[jj], pch = 20)
          points(antherPoll$anthx[jj], antherPoll$anthy[jj], pch = 20, col = 'red')
         #Sys.sleep(0.1)
     
     }
     
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
     
     
     plot(antherPoll$anthx.1.abs, antherPoll$anthy.1.abs)
     f = fitSS(na.omit(cbind(antherPoll$anthx.1.abs, antherPoll$anthy.1.abs)))
     
     
     efit <- fit.ellipse(na.omit(antherPoll$anthx.1.abs), na.omit(antherPoll$anthy.1.abs))
     e <- get.ellipse(efit)
     plot(e, col = 'red', type = 'b')
     points(antherPoll$anthx.1.abs, antherPoll$anthy.1.abs, asp = 1)
     
     points(efit$center[1], efit$center[2])
     draw.circle(x = f$par[1], y = f$par[2] + 0.001,radius = f$par[3])



}




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
]