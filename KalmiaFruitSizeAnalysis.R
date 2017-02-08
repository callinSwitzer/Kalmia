# Kalmia analyze fruit size at end of season
# Using fruit sizes calculated from image segmentation in Python with opencv
#  Callin Switzer
# 20 October 2016
# 8 Feb 2017 Update: Conducted Statistical Modeling with LMER and GLMER



ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "gsheet", 'lme4', 'plyr')
ipak(packages)

theme_set(theme_classic())


URL <- 'https://docs.google.com/spreadsheets/d/1FyN3Pnks3CP3c8MA26qZt70OKK5a0dn9TufnhEUkTMI/edit?usp=sharing'

# read in data
kfrt <- gsheet2tbl(URL)


# clean and process data
kfrt <- kfrt[kfrt$dia_mm != 20.0, ]

nrow(kfrt)

kfrt$trt <- sapply(X = 1:nrow(kfrt), FUN = function(x) strsplit(kfrt$plantNum[x], split = "__")[[1]][2])
kfrt$accessNum <- sapply(X = 1:nrow(kfrt), FUN = function(x) strsplit(kfrt$plantNum[x], split = "__")[[1]][1])

kfrt$trt <- mapvalues(kfrt$trt, c(1,2,3,4,5), c("Bagged", "Bagged & Selfed", "Unbagged", 
                                                "Unbagged & Outcrossed", "Unbagged_2"))



# count number of fruits
counts = as.data.frame(table(kfrt$plantNum))
counts$trt = sapply(X = 1:nrow(counts), FUN = function(x) strsplit(as.character(counts$Var1[x]), split = "__")[[1]][2])

counts$trt <- mapvalues(counts$trt, c(1,2,3,4,5), c("Bagged", "Bagged & Selfed", "Unbagged", 
                                                "Unbagged & Outcrossed", "Unbagged_2"))


counts$accessNum = sapply(X = 1:nrow(counts), FUN = function(x) strsplit(as.character(counts$Var1[x]), split = "__")[[1]][1])


# add in the trts and accession numbers that had counts of 0
kalNotes <- gsheet2tbl('https://docs.google.com/spreadsheets/d/1IidOr3GP77RCt5vRUe60FHIpy4emx7NwE9JQGqmUzjI/edit?usp=sharing')

unique(kfrt$plantNum)
# get the accession numbers I should have
accNumsHave <- unique(kalNotes$Plant.Number)
accNumsHave <- accNumsHave[!(accNumsHave %in% c("", "Plant Number"))]

#change formatting to match above
accNumsHave <- gsub("#", "", accNumsHave)
accNumsHave <- gsub("\ ", "_", accNumsHave)
accNumsHave <- gsub("\\-", "_", accNumsHave)
accNumsHave <- gsub("\\*", "_", accNumsHave)
accNumsHave <-toupper(accNumsHave)

shouldHave <- as.data.frame(t(sapply(accNumsHave, function(x) as.character(paste(x, c(1:4), sep = "__")))))
row.names(shouldHave) <- NULL
shouldHave1 <- c(as.character(shouldHave[, 1]), 
                 as.character(shouldHave[, 2]), 
                 as.character(shouldHave[, 3]), 
                 as.character(shouldHave[, 4]))

# find missing ones -- this is the ones that
# are in the daily datasheet that aren't in the fruit measurements
missingTrts <- setdiff(shouldHave1,unique(kfrt$plantNum[kfrt$trt != "5"]))

# this should be 0 -- the ones from the fruit measurements that aren't in the daily datasheet
setdiff(unique(kfrt$plantNum[kfrt$trt != "5"]), shouldHave1)

# here's the list of plants that had their bags/tags go missing during the experiment (meaning that
# I was unable to collect fruits)
# note: "677_66_MASS__1" was the only plant that was missing a tag during the final collection
# that wasn't missing in one of my previous checks. 
NAList <- c("1129_74_E__4", "1129_74_C__4", "39_60_A__3", "685_70_A__2", "677_66_MASS__1")

# note: this ignores trt#5 which is the same as #3
ZeroFruitPlants <- missingTrts[!(missingTrts %in% NAList)]

zfdf <- data.frame(Var1 = ZeroFruitPlants, 
                   Freq = 0, 
                   trt = sapply(X = 1:length(ZeroFruitPlants), 
                                FUN = function(x) strsplit(as.character(ZeroFruitPlants[x]), 
                                                           split = "__")[[1]][2]),
                   accessNum = sapply(X = 1:length(ZeroFruitPlants), 
                                            FUN = function(x) strsplit(as.character(ZeroFruitPlants[x]), 
                                                                       split = "__")[[1]][1])
                   )

zfdf$trt <- mapvalues(zfdf$trt, c(1,2,3,4,5), c("Bagged", "Bagged & Selfed", "Unbagged", 
                                                "Unbagged & Outcrossed", "Unbagged_2"))

countFin <- rbind(counts, zfdf)


# final fruit counts for the fruits collected at the end of the experiment
countFin <- countFin[order(countFin$accessNum, countFin$trt, decreasing = FALSE), ]


# visualize counts of fruits
ggplot(countFin, aes(x = trt , y = Freq, fill = trt)) + 
     geom_boxplot() +
    # geom_violin()+
     
     labs(y = "Number of fruits", x = "Treatment") + 
     scale_fill_brewer(name = "Treatment", palette = "Set1")

summary(aov(countFin$Freq ~ countFin$trt))
saveDir <- "/Users/callinswitzer/Dropbox/ExperSummer2016/Kalmia/"
ggsave(paste0(saveDir, "KalmiaFruitNumber.pdf"), width = 10, height = 8)

# visualize counts of fruits (just for 4 trts)
ggplot(countFin[!(countFin$trt %in% 'Unbagged_2'), ], aes(x = trt , y = Freq, fill = trt)) + 
     geom_boxplot() +
     labs(y = "Number of fruits", x = "Treatment") +
     theme(axis.text.x = element_text(angle = 35, hjust = 1), 
           legend.position = 'none') + 
     scale_fill_brewer(name = "Treatment", palette = "Set1")

ggsave(paste0(saveDir, "KalmiaFruitNumber_trt1_4Only.pdf"), width = 5, height = 4)

# visualize fruit size
sizeDF_mean <- as.data.frame(tapply(kfrt$dia_mm, INDEX = kfrt$plantNum, mean))
colnames(sizeDF_mean) = "meanFrtSz"
sizeDF_mean$trt = sapply(X = 1:length(sizeDF_mean$meanFrtSz), 
                          FUN = function(x) strsplit(as.character(row.names(sizeDF_mean)[x]), 
                                                     split = "__")[[1]][2])

sizeDF_mean$trt <- mapvalues(sizeDF_mean$trt, c(1,2,3,4,5), c("Bagged", "Bagged & Selfed", "Unbagged", 
                                                "Unbagged & Outcrossed", "Unbagged_2"))
sizeDF_mean$accessNum = sapply(X = 1:length(sizeDF_mean$meanFrtSz), 
                   FUN = function(x) strsplit(as.character(row.names(sizeDF_mean)[x]), 
                                              split = "__")[[1]][1])


ggplot(sizeDF_mean, aes(x = trt, y = meanFrtSz, fill = trt)) + 
     geom_boxplot(alpha = 0.5) + 
#      stat_summary(fun.y=mean, geom="line", aes(group = accessNum, color = accessNum))  + 
#      stat_summary(fun.y=mean, geom="point", aes(group = accessNum, color = accessNum)) + 
     geom_point(aes(color = trt))


ggplot(sizeDF_mean, aes(x = trt, y = meanFrtSz, fill = trt)) + 
     geom_boxplot() + 
     labs(y = "Mean Fruit Diameter (mm)", x = "Treatment") + 
     scale_fill_brewer(name = "Treatment", palette = "Set1")
ggsave(paste0(saveDir, "KalmiaFruitDiameter.pdf"), width = 10, height = 8)


ggplot(sizeDF_mean[!(sizeDF_mean$trt %in% 'Unbagged_2'), ], 
       aes(x = trt, y = meanFrtSz, fill = trt)) + 
     geom_boxplot() + 
     labs(y = "Mean Fruit Diameter (mm)", x = "Treatment") + 
     scale_fill_brewer(name = "Treatment", palette = "Set1") + 
     theme(axis.text.x = element_text(angle = 35, hjust = 1), 
          legend.position = 'none') 
ggsave(paste0(saveDir, "KalmiaFruitDiameter_trt14Only.pdf"), width = 5, height = 4)

# compare fruit size with number of seeds
sizeSeed <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1pm7-1HD5fnhyMmV-StPh_geE44MJxT5thC_AlfRQ0U4/edit?usp=sharing")

ggplot(sizeSeed, aes(x = Dia..mm., y = NumSeeds)) + 
     geom_point() + 
     stat_smooth(method = 'lm', formula = y ~ exp(x), se = F) + 
     labs(x = 'Fruit Diameter (mm)', y = 'Num Seeds in 1/5 of Fruit')

ggsave(paste0(saveDir, "KalmiaFruitSeeds.pdf"), width = 5, height = 4)


# on log scale
ggplot(sizeSeed, aes(x = Dia..mm., y = NumSeeds)) + 
     geom_point() + 
     stat_smooth(method = 'lm', formula = y ~ x, se = F) + 
     scale_y_continuous(trans="log") + 
     labs(y = "log(e) number of seeds")
