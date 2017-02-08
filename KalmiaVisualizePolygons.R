## import polygons of pollen area from matlab
## align them and plot them 

# read file names
library(gsheet)

# KalmiaAndDeadBees
kds <- gsheet2tbl("https://docs.google.com/spreadsheets/d/17YTmbErH53yR0AYLfwWCmTocqcrqLQFqfj2e8kT-pig/edit?usp=sharing")

unique(kds$FlowerAccNum)

# directory where the polygon files are found
setwd('/Users/callinswitzer/Dropbox/KalmiaProject/KalmiaExamples/KalmiaManualTrig_VidsToProcess/')


# Here start


polyFiles <- 



p1 <- read.table("~/Desktop/20160616_154351_polygon.txt", sep = ',', header = TRUE)
p2 <- read.table("~/Desktop/20160615_163530_polygon.txt", sep = ',', header = TRUE)
p3 <- read.table("~/Desktop/20160615_162452_polygon.txt", sep = ',', header = TRUE)
p4 <- read.table("~/Desktop/20160616_155529_polygon.txt", sep = ',', header = TRUE)


p1$pg <- "one"
p2$pg <- "two"
p3$pg <- "three"
p4$pg <- "four"

pgs <- rbind(p1, p2, p3, p4)
plot(x = p1$xvals, y = 600-p1$yvals, type = 'n', xlim = c(-200, 600))
polygon(x = p1$xvals[1:(nrow(p1) - 4) ], y = 600-p1$yvals[1:(nrow(p1) - 4) ], col = rgb(0,0,1, 0.3), border = NA)

polygon(x = p2$xvals[1:(nrow(p2) - 4) ]-35, y = 600-p2$yvals[1:(nrow(p2) - 4) ], col = rgb(0,0,1, 0.3), border = NA)


polygon(x = p3$xvals[1:(nrow(p3) - 4) ]-35, y = 600-p3$yvals[1:(nrow(p3) - 4) ] + 30, col = rgb(0,0,1, 0.3), border = NA)

polygon(x = p4$xvals[1:(nrow(p4) - 4) ]-35, y = 600-p4$yvals[1:(nrow(p4) - 4) ] + 30, col = rgb(0,0,1, 0.3), border = NA)



polygon(x = 440-p1$xvals[1:(nrow(p1) - 4) ], y = 600-p1$yvals[1:(nrow(p1) - 4) ], col = rgb(0,0,1, 0.3), border = NA)

polygon(x = 500-p2$xvals[1:(nrow(p2) - 4) ]-35, y = 600-p2$yvals[1:(nrow(p2) - 4) ], col = rgb(0,0,1, 0.3), border = NA)


polygon(x = 500-p3$xvals[1:(nrow(p3) - 4) ]-35, y = 600-p3$yvals[1:(nrow(p3) - 4) ] + 30, col = rgb(0,0,1, 0.3), border = NA)

polygon(x = 500-p4$xvals[1:(nrow(p4) - 4) ]-35, y = 600-p4$yvals[1:(nrow(p4) - 4) ] + 30, col = rgb(0,0,1, 0.3), border = NA)

# rotate polygons
#read in X and Y as vectors
X = p1$xvals[1:(nrow(p1) - 4) ]
Y = 600-p1$yvals[1:(nrow(p1) - 4) ]

X = p1$xvals
Y = 600-p1$yvals
M <- cbind(X,Y)
#plot data
plot(M[,1],M[,2],xlim=c(0,1200),ylim=c(0,1200), type = 'l')
#calculate rotation angle
x1 <- M[(nrow(p1) - 2),1 ]
y1 <- M[(nrow(p1) - 2),2 ]
x2 <- M[(nrow(p1)),1 ]
y2 <- M[(nrow(p1)),2 ]

points(c(x1, x2), c(y1, y2))


alpha <- -atan((y1-y2)/(x1-x2)) + pi/2
#rotation matrix
rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
#shift, rotate, shift back
M2 <- t(rotm %*% (
     t(M)-c(M[1,1],M[1,2]))+c(M[1,1],M[1,2]))

x3 <-  M2[(nrow(p1) - 2),1 ]
y3 <- M2[(nrow(p1) - 2),2 ]
x4 <-M2[(nrow(p1)),1 ]
y4 <- M2[(nrow(p1)),2 ]

#plot
plot(M2[,1],M2[,2],xlim=c(0,1200),ylim=c(0,1200), type = 'l')
points(x = c(x3, x4), y= c(y3, y4))


# move
M2[,1] <- M2[, 1] - x4
M2[,2 ] <- M2[,2 ] - y4
par(pty="s")
plot(M2[,1],M2[,2],xlim=c(-600,600),ylim=c(0,1200), type = 'l')

# scale 
M2 <- M2 * 1.5
lines(M2[,1],M2[,2],xlim=c(0,1200),ylim=c(0,1200), type = 'l')



plot(M2[,1],M2[,2],xlim=c(0,1200),ylim=c(0,1200), type = 'l')



#### 
xx <- read.csv('~/Downloads/graphx.csv', header = F)
yy <- read.csv('~/Downloads/graphy.csv', header = F)
M <- cbind(xx[,1], yy[,1])

#read in X and Y as vectors
#plot data
plot(M[,1],M[,2],xlim=c(0,1200),ylim=c(0,1200))
#calculate rotation angle

x1 <-  M[1,1]
y1 <- M[1,2]
x2 <- tail(M,1)[,1]
y2 <- tail(M,1)[,2]


alpha <- -atan((y1-y2)/(x1-x2)) 
#alpha <- -atan((M[1,2]-tail(M,1)[,2])/(M[1,1]-tail(M,1)[,1]))
#rotation matrix
rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
#shift, rotate, shift back
M2 <- t(rotm %*% (
     t(M)-c(M[1,1],M[1,2])
)+c(M[1,1],M[1,2]))
#plot
points(M2[,1],M2[,2],xlim=c(0,1200),ylim=c(0,1200))
lines(y = c(M2[1,2], tail(M2,1)[,2]), x = c(M2[1,1],tail(M2,1)[,1]), pch = 20)
