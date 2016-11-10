## import polygons of pollen area from matlab
## align them and plot them 

# read file names
 library(gsheet)

#file.choose()
setwd("/Users/callinswitzer/Dropbox/KalmiaProject/KalmiaExamples/KalmiaManualTrig_VidsToProcess/")

# get files with polygon coordinates
polyFiles <- dir()[grep(dir(), pattern = "_polygon.txt") ]

## read in files, and make into a long dataframe 
## remove last 9 points which contain metadata

for(ii in 1:length(polyFiles)){
     tmp <- read.csv(polyFiles[ii])
     #read in X and Y as vectors -- these are the four points on the flower for reference
     # rotate polygons
     X = tmp$xvals
     Y = 600-tmp$yvals
     M <- cbind(X,Y)
     #plot data
     # plot(M[,1],M[,2],xlim=c(-600,600),ylim=c(-100,1200), type = 'l')
     #calculate rotation angle
     x1 <- M[(length(X) - 7),1 ]
     y1 <- M[(length(Y) - 7), 2 ]
     x2 <- M[(length(X)-5),1 ]
     y2 <- M[(length(Y)-5),2 ]
     
     # points(c(x1, x2), c(y1, y2))
     
     
     alpha <- -atan((y1-y2)/(x1-x2)) + pi/2
     #rotation matrix
     rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
     
     
     #shift points, so that turning point is (0,0)
     M2.1 <- t(t(M)-c(x1,y1))
     # polygon(M2.1,col="blue")
     # abline(h = 0, v =0)
     
     #rotate
     M2.2 <- t(rotm %*% (t(M2.1)))
     # polygon(M2.2,col="green")
     x12 <- M2.2[(length(X) - 7),1 ]
     y12 <- M2.2[(length(Y) - 7), 2 ]
     x22 <- M2.2[(length(X)-5),1 ]
     y22 <- M2.2[(length(Y)-5),2 ]
     
     # points(c(x12, x22), c(y12, y22))
     # plot(tail(M2.2,  n = 9 ))
     
     ## make sure that the points aren't upside down (rotate 180 degrees, if necessary)
     if(y22 > y12){
          alpha <-  pi
          #rotation matrix
          rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
          M2.2 <- t(rotm %*% (t(M2.2)))
     }
     
     # scale so that length between center of flower and pocket is 1
     x13 <- M2.2[(length(X) - 8),1 ]
     y13 <- M2.2[(length(Y) - 8), 2 ]
     x23 <- M2.2[(length(X)-7),1 ]
     y23 <- M2.2[(length(Y)-7),2 ]
     
     # points(c(x13, x23), c(y13, y23), col = 'red', pch = 20)
     
     # find distance 
     dst <- sqrt((x13-x23)^2 + (y13 -y23)^2)
     
     M2.3 <- M2.2 / dst
     # plot(M2.3, type ='n')
     # polygon(M2.3, type ='n')
     
#      if(ii == 1) {
#           plot(M2.3[(nrow(M2.3) - 8):nrow(M2.3),])
#           text(M2.3[(nrow(M2.3) - 8):nrow(M2.3),], labels = 1:9, adj =  -0.5)
#      }
#      else {
#           points(M2.3[(nrow(M2.3) - 8):nrow(M2.3),])
#           text(M2.3[(nrow(M2.3) - 8):nrow(M2.3),], labels = 1:9, adj =  -0.5)
#      }

     
     tmp <- as.data.frame(M2.3[1:(nrow(M2.3) - 9),])
     tmp$rec = substr(polyFiles[ii], start = 1, stop = 15)
     if(ii == 1) tmp2 <- tmp
     else tmp2 <- rbind(tmp, tmp2)
     
     
     if (ii == 1){
          plot(tmp[, 1], tmp[,2], type = 'n', xlim = c(-3,3), ylim = c(-1, 5))
          polygon(tmp[, 1], tmp[,2], col = rgb(0,0,1,0.05), border = NA)

     } 
     else {
          ## To flip horizontaly ,randomly
          tmp[,1] <- sample(x = c(1, -1), size = 1) * tmp[,1]
          polygon(tmp[, 1], tmp[,2], col = rgb(0,0,1,0.05), border = NA)
     }
     Sys.sleep(0.1)
}

# save polygons as images

for(ii in 1:length(polyFiles)){
     tmp <- read.csv(polyFiles[ii])
     #read in X and Y as vectors -- these are the four points on the flower for reference
     # rotate polygons
     X = tmp$xvals
     Y = 600-tmp$yvals
     M <- cbind(X,Y)
     #plot data
     # plot(M[,1],M[,2],xlim=c(-600,600),ylim=c(-100,1200), type = 'l')
     #calculate rotation angle
     x1 <- M[(length(X) - 7),1 ]
     y1 <- M[(length(Y) - 7), 2 ]
     x2 <- M[(length(X)-5),1 ]
     y2 <- M[(length(Y)-5),2 ]
     
     # points(c(x1, x2), c(y1, y2))
     
     
     alpha <- -atan((y1-y2)/(x1-x2)) + pi/2
     #rotation matrix
     rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
     
     
     #shift points, so that turning point is (0,0)
     M2.1 <- t(t(M)-c(x1,y1))
     # polygon(M2.1,col="blue")
     # abline(h = 0, v =0)
     
     #rotate
     M2.2 <- t(rotm %*% (t(M2.1)))
     # polygon(M2.2,col="green")
     x12 <- M2.2[(length(X) - 7),1 ]
     y12 <- M2.2[(length(Y) - 7), 2 ]
     x22 <- M2.2[(length(X)-5),1 ]
     y22 <- M2.2[(length(Y)-5),2 ]
     
     # points(c(x12, x22), c(y12, y22))
     # plot(tail(M2.2,  n = 9 ))
     
     ## make sure that the points aren't upside down (rotate 180 degrees, if necessary)
     if(y22 > y12){
          alpha <-  pi
          #rotation matrix
          rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
          M2.2 <- t(rotm %*% (t(M2.2)))
     }
     
     # scale so that length between center of flower and pocket is 1
     x13 <- M2.2[(length(X) - 8),1 ]
     y13 <- M2.2[(length(Y) - 8), 2 ]
     x23 <- M2.2[(length(X)-7),1 ]
     y23 <- M2.2[(length(Y)-7),2 ]
     
     # points(c(x13, x23), c(y13, y23), col = 'red', pch = 20)
     
     # find distance 
     dst <- sqrt((x13-x23)^2 + (y13 -y23)^2)
     
     M2.3 <- M2.2 / dst
     # plot(M2.3, type ='n')
     # polygon(M2.3, type ='n')
     
     #      if(ii == 1) {
     #           plot(M2.3[(nrow(M2.3) - 8):nrow(M2.3),])
     #           text(M2.3[(nrow(M2.3) - 8):nrow(M2.3),], labels = 1:9, adj =  -0.5)
     #      }
     #      else {
     #           points(M2.3[(nrow(M2.3) - 8):nrow(M2.3),])
     #           text(M2.3[(nrow(M2.3) - 8):nrow(M2.3),], labels = 1:9, adj =  -0.5)
     #      }
     
     
     tmp <- as.data.frame(M2.3[1:(nrow(M2.3) - 9),])
     tmp$rec = substr(polyFiles[ii], start = 1, stop = 15)
     if(ii == 1) tmp2 <- tmp
     else tmp2 <- rbind(tmp, tmp2)
     
     
     # plot polygon and save
     # save as high resolution
     png( paste0(tmp$rec, 'TEST.png'), width = 2000, height = 2000)
     par(mai=c(0,0,0,0))
     par( mgp = c(0, -1.4, 0)  )
     par(tcl = 0.5)

     plot(tmp[, 1], tmp[,2], type = 'n', xlim = c(-4,4), ylim = c(-1, 7), 
          asp=1, ann=FALSE, bty = 'n')
     polygon(tmp[, 1], tmp[,2], col = rgb(0,0,1,1), border = NA)
     #lines(x = c(.1,.1, -.1, -.1, .1, -1, 1), y = c(0, -.5, 0, -.5, 0,0, 0))
     abline(h = 0, v = c(0, 1, -1))
     dev.off()
     
     
}


# for(jj in 1:30){
#      # bootstrap heatmap
#      rec_boot <- sample(unique(tmp2$rec), replace = TRUE)
#      
#      # plot all polygons
#      for(ii in 1:length(rec_boot)){
#           tmp <- tmp2[tmp2$rec == rec_boot[ii], ]
#           if (jj == 1){
#                plot(tmp[, 1], tmp[,2], type = 'n', xlim = c(-3,3), ylim = c(-1, 5))
#                polygon(tmp[, 1], tmp[,2], col = rgb(0,0,1,0.002), border = NA)
#                
#           } 
#           else {
#                ## To flip horizontaly ,randomly
#                tmp[,1] <- sample(x = c(1, -1), size = 1) * tmp[,1]
#                polygon(tmp[, 1], tmp[,2], col = rgb(0,0,1,0.002), border = NA)
#           }
#      }
#      
# }





