# Setup
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "plyr", "car", "lme4", "plotrix")

ipak(packages); rm(packages)

# example from here: http://cfsv.synechism.org/c2/sec23.pdf


# example data
fps = 100
tme = 1:500  / fps # seconds
position = data.frame(xx = 2 * cos(tme), yy = sin(tme))
dim(position)
quartz(title = "position")
par(mfrow =c(2,2))
plot(position)
text(x = 2, y = 0, labels = "Start", adj = c(0.5, 1.5))

plot(position$xx)
plot(position$yy)


# velocity
velocity = data.frame(rbind(c(NA, NA), apply(position, MARGIN = 2, FUN = diff))) * fps
quartz("velocity")
par(mfrow =c(2,2))
speed = sqrt(velocity$xx^2 + velocity$yy^2)
plot(velocity$xx, x = tme) # calculated
lines(y = -2*sin(tme), x = tme, col = 'red') # theoretical
plot(velocity$yy, x = tme)
lines(y = cos(tme), x = tme, col = 'red') # theoretical

plot(speed, x = tme)
lines(y = sqrt(3 * sin(tme)^2 + 1), x = tme, col = 'red') # theoretical 

# acceleration
acceleration = data.frame(rbind(c(NA, NA), apply(velocity, MARGIN = 2, FUN = diff))) * fps
quartz("acceleration")
par(mfrow =c(2,2))
plot(acceleration$xx, x = tme) # calculated
lines(y = -2*cos(tme), x = tme, col = 'red') # theoretical
plot(acceleration$yy, x = tme)
lines(y = -sin(tme), x = tme, col = 'red') # theoretical


# unit tangent vector
T_t = velocity / speed
quartz("UnitTangentVector")
par(mfrow = c(2,2))
plot(T_t$xx, x = tme)  # calculated
lines(sqrt(2/(5-3*cos(2 * tme))) * -2*sin(tme), x = tme, col = 'red') # theoretical
plot(T_t$yy, x = tme)
lines(sqrt(2/(5-3*cos(2 * tme))) * cos(tme), x = tme, col = 'red') # theoretical

DT = data.frame(rbind(c(NA, NA), apply(T_t, MARGIN = 2, FUN = diff))) * fps
plot(DT$xx, x = tme) # calculated
lines(y = sqrt(2/(5-3*cos(2 * tme))) * 
           -2*cos(tme) - (3*sqrt(2) * sin(2*tme) 
                          / (5 - 3*cos(2*tme))^(3/2)) 
      * -2*sin(tme), x = tme, col= 'red') # theoretical

plot(DT$yy, x = tme)
lines(y = sqrt(2/(5-3*cos(2 * tme))) * 
           -sin(tme) - (3*sqrt(2) * sin(2*tme) 
                          / (5 - 3*cos(2*tme))^(3/2)) 
      * cos(tme), x = tme, col= 'red') # theoretical

# check time point 0.79 (pi / 4)
position[tme == 0.79, ]
sqrt(2) #theory
1/sqrt(2)

velocity[tme == 0.79, ]

speed[tme == 0.79]
sqrt(5/2)#theory

T_t[tme == 0.79, ]
1/sqrt(5) * -2#theory
1/sqrt(5) * 1

DT1 = DT[tme == 0.79, ]
aa = -1 / (5 * sqrt(5)) * 4 #theory
bb = -1 / (5 * sqrt(5)) * 8

NormDT = sqrt(DT$xx^2 + DT$yy^2)
4/5 # theory

Curvature = NormDT / speed
Curvature[tme == 0.79]
(4/5) / sqrt(5/2) # theory


# compute a_N (normal acceleration) and a_T (tangential acceleration)
# we know a

# a_T = ds/dt
a_T = dsdt = c(NA, diff(speed) * fps)
dsdt[tme == 0.79] # tangential acceleration calculation
1/sqrt(2) * (5 - 3 * cos(2 * pi / 4))^(-0.5)* (3* sin(2 * pi/4)) # theory tangential accel
3 / sqrt(10) # theory

# a_N = speed^2 * curvature
a_N = speed^2 * Curvature
a_N[tme == 0.79] # calculated
4/sqrt(10) # theoretical -- great!

# check a = a



acceleration[tme == 0.79, ]
-sqrt(2) * 1/sqrt(5) * -2
-1/(sqrt(2)) * 1/sqrt(5) * 1

sqrt((-sqrt(2) * 1/sqrt(5) * -2)^2 + (-1/(sqrt(2)) * 1/sqrt(5) * 1)^2)

3/sqrt(10)

# dot product version of calculation of tangential acceleration
#a_T = acceleration dot T_t
c(-sqrt(2), -1/sqrt(2)) %*% (1/sqrt(5) * c(-2, 1)) # theory
a_T_ver2 <- sapply(X = 1:nrow(acceleration), FUN  = function(x) unlist(acceleration[x, ]) %*% unlist(T_t[x, ]))
a_T_ver2[tme == 0.79] # calculation close
a_T[tme == 0.79] # other calculation
3/sqrt(10) #theory


# dot product version of calculating normal acceleration
# a_N = acceleration dot N_t
# where N_t = DT  / ||DT|| where, ||DT|| is the norm of DT
N_t = data.frame(t(sapply(1:nrow(DT), FUN = function(x) unlist(DT[x, ] / NormDT[x]))))
a_N_ver2 <- sapply(X = 1:nrow(acceleration), FUN  = function(x) unlist(acceleration[x, ]) %*% unlist(N_t[x,]))
a_N_ver2[tme == 0.79] # calculated with one method
a_N[tme == 0.79] # calculated
4/sqrt(10) # theoretical -- close

# check total accel by adding normal and tangential accelerations
# a_total = a_T * T_t + a_N * N_t
a_total = as.data.frame(t(sapply(X = 1:nrow(acceleration), FUN  = function(x) a_T[x] * T_t[x, ] + a_N[x] * N_t[x,] )))
plot(a_total)
plot(acceleration)
quartz("atotal")
par(mfrow = c(2,2))
plot(unlist(a_total$xx))
plot(unlist(a_total$yy))
plot(acceleration$xx)
plot(acceleration$yy)

plot(a_N)
plot(a_T)


plot(a_T, x = tme)
plot(a_N, x = tme)

T_t[tme == 0.79, ]
1/sqrt(5) * c(-2, 1)

a_T[tme == 0.79]

# calculate magnitude of acceleration, using two methods
# 1. Normal and tangential acceleration
a_mag1 = sqrt(a_T^2 + a_N^2)
plot(a_mag1)

amag2 = sqrt(acceleration$xx^2 + acceleration$yy^2)
plot(amag2)

round(amag2, 1) == round(a_mag1, 1)

df1 = data.frame(amag2, a_mag1, abs(amag2 - a_mag1))
df1

# Questions for Stacey:
# Do people report magnitude of acceleration by combining X and Y
# or do they report tangential and normal acceleration separately?



# BUtterworth filtering
library(signal)
n <- 100
x <- 1:n
y <- ifelse(0.3*n < x & x < 0.7*n, 1, 0)
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.7, 0))
plot(x, y, type='o', pch=20, ylim=c(-0.1, 1.1))
W <- 0.5 # 1/2 the nyquist frequency
b1 <- butter(1, W)
y1 <- filtfilt(b1, y)
points(x, y1, type='o', pch=20, col='red')

b2 <- butter(2, W)
y2 <- filtfilt(b2, y)
points(x, y2, type='o', pch=20, col='blue')

b3 <- butter(5, W)
y3 <- filtfilt(b3, y)
points(x, y3, type='o', pch=20, col='forestgreen')


sg <- sgolayfilt(y, p = 1, n = 41) # Savitzky-Golay filter
plot(y, type="b", col = 'red', pch = 20)
lines(sg, pch = 20, type = 'l') # smoothed SG data

legend("topright", lwd=2, pch=20, 
       col=c("black", "red", "blue", "forestgreen"),
       legend=c("Signal", "Order 1", "Order 2", "Order 3"))
