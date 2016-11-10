# Callin Switzer
# 7/28/2016
# Analyze and visualize the kalmia from the Arboretum crossing 
# experiments


ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "scales", "multcomp", "contrast", "plyr", "pwr", "gsheet", 'lme4')
ipak(packages)


url <- "https://docs.google.com/spreadsheets/d/1IidOr3GP77RCt5vRUe60FHIpy4emx7NwE9JQGqmUzjI/edit?usp=sharing"
kalm <- gsheet2tbl(url)

# remove extra rows
kalm <- kalm[!(kalm$Location %in% c("Location", "")),]

# convert to date
kalm$Date <- as.POSIXlt(kalm$Date, format = "%d %b %Y", tz = "EST")

colnames(kalm)


library(reshape2)

kalm.long1<-reshape(kalm[as.character(kalm$Date) == "2016-07-04", ], varying=c("buds_flowers.1", "buds_flowers.2", "buds_flowers.3", "buds_flowers.4", "flowers.1", "flowers.2", "flowers.3", "flowers.4", "fruits.1", "fruits.2", "fruits.3", "fruits.4"), direction="long", idvar="Plant.Number", sep=".", timevar = "treatment")


ggplot(kalm, aes(x = Date, y = as.numeric(buds_flowers.1), 
                 color = Plant.Number)) + 
     geom_line(size =2, alpha = 0.5) + 
     theme_bw()

ggplot(kalm[kalm$Plant.Number == "150-58*A", ], aes(x = Date, y = as.numeric(buds_flowers.1), color = Plant.Number)) + 
     geom_line(size =2, alpha = 0.5) + 
     theme_bw()
unique(kalm$Plant.Number)

ggplot(kalm, aes(x = Date, y = as.numeric(fruits.4), color = Plant.Number)) + 
     geom_line(size =2, alpha = 0.5) + 
     theme_bw()

ggplot(kalm[kalm$Plant.Number == "1129-74*B", ], aes(x = Date, y = as.numeric(fruits.4), color = Plant.Number)) + 
     geom_line(size =2, alpha = 0.5) + 
     theme_bw()


# combine into one long data frame
ii <- 5
kl <- cbind(kalm[, 1:4], count = kalm[,ii])
kl$measurement <- strsplit(as.character(colnames(kalm)[ii]), split = "\\.")[[1]][1]
kl$treatment <- strsplit(as.character(colnames(kalm)[ii]), split = "\\.")[[1]][2]

kl$treatment <- mapvalues(kl$treatment, c(1,2,3,4), c("Bagged", "Bagged & Selfed", "Unbagged", 
                                                "Unbagged & Outcrossed"))


colnames(kalm)[1:25]
nums <- c(6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19)
for(ii in nums){
     foo <- cbind(kalm[, 1:4], count = kalm[,ii])
     foo$measurement <- strsplit(as.character(colnames(kalm)[ii]), split = "\\.")[[1]][1]
     foo$treatment <- strsplit(as.character(colnames(kalm)[ii]), split = "\\.")[[1]][2] 
     
     kl <- rbind(kl, foo)
}

kl$count <- as.numeric(kl$count)
kl$treatment <- mapvalues(kl$treatment, c(1,2,3,4), c("Bagged", "Bagged & Selfed", "Unbagged", 
                                                      "Unbagged & Outcrossed"))




ggplot(kl[kl$measurement == "flowers",], aes(x = Date, y = count, color = treatment)) + 
     geom_smooth(se = F, span = 0.1) + 
     theme_bw() + 
     labs(x = "Date", y = "Number of Open Flowers") + 
     scale_color_brewer(name = "Treatment", palette = "Set1")
saveDir <- "/Users/callinswitzer/Dropbox/ExperSummer2016/Kalmia/"
ggsave(paste0(saveDir, "KalmiaOpenFlowers.pdf"), width = 6, height = 4)

ggplot(kl[kl$measurement == "fruits",], aes(x = Date, y = count, color = treatment)) + 
     #geom_point() + 
     geom_smooth(se = F, span = 0.2, size = 2) + 
     theme_bw() + 
     labs(x = "Date", y = "Number of Fruits") + 
     scale_color_brewer(name = "Treatment", palette = "Set1")
ggsave(paste0(saveDir, "KalmiaNumFruits.pdf"), width = 6, height = 4)


ggplot(kl[kl$measurement == "buds_flowers",], aes(x = Date, y = count, color = treatment)) + 
     geom_smooth(se = F, span = 0.1, size = 2) + 
     theme_bw() + 
     labs(y = "Num Buds + Flowers") + 
     scale_color_brewer(name = "Treatment", palette = "Set1")
ggsave(paste0(saveDir, "KalmiaBudsFlowers.pdf"), width = 6, height = 4)

unique(kl$measurement)

ggplot(kl[kl$measurement == "flowers_w_flipped_anth",], aes(x = Date, y = count, color = treatment)) + 
     geom_smooth(se = F, span = 0.1, size = 2) + 
     theme_bw() + 
     labs(y = "Number of flowers w/ at least 1 anther flipped") + 
     scale_color_brewer(name = "Treatment", palette = "Set1")
ggsave(paste0(saveDir, "KalmiaFlippedAnthers.pdf"), width = 6, height = 4)

