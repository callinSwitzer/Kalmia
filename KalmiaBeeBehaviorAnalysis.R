# Callin Switzer
# 9 Nov 2016
## Analyze the behavior of pollinators
## from videos that I manually classified
## added to github

set.seed(12345)
ipak <- function(pkg){
     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
     if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
     sapply(pkg, require, character.only = TRUE)
}

packages <- c("ggplot2", "gsheet", 'lme4', 'plyr', 'nnet', 'viridis')
ipak(packages)

# set ggplot theme
theme_set(theme_classic())


# gsheet called Kalmia Bee Behavior classification
URL = 'https://docs.google.com/spreadsheets/d/1gxAmCcwJ9zsnYxjXe2hnTowPUxTjadbBx-eXJYS4_Js/edit?usp=sharing'

beh = gsheet2tbl(URL)

# discard some messed up data (labeled as discard)
discRows = grep(pattern = 'discard', 
                x = tolower(apply(beh, 
                         MARGIN = 1, 
                         FUN = paste0, collapse = "_")))
beh = beh[-discRows, ]

# remove lab views with captive bees
beh_field <-beh[beh$vid.type == 'field',]


# visualize diversity of pollinators
colnames(beh_field)
unique(beh_field$pollinatorClass)

# set order of levels
beh_field <- within(beh_field, 
                    pollinatorClass <- factor(pollinatorClass, 
                            levels=names(sort(table(pollinatorClass), 
                                              decreasing=TRUE))))
ggplot(beh_field, aes(x = pollinatorClass)) + 
     geom_bar() + 
     labs (x = "Insect visitor", y = "Frequency") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
           axis.ticks.x=element_blank()) + 
     scale_x_discrete(labels=c("Bombus spp." = expression(
          paste(italic("Bombus spp."))), 
          "Xylocopa virginica" = expression(
               paste(italic("Xylocopa virginica"))),
          "Apis mellifera" = expression(
               paste(italic("Apis mellifera")))))
ggsave(filename = "KalmiaFigures/VisitorProps.pdf", width = 5, height = 4)

# make proportion to show figure

propTab <- data.frame(prop.table(table(beh_field$pollinatorClass)))
names(propTab) = c("Insect", "Proportion")
propTab

# resample proportion table to get error bars
resTab <- function(o){
     prop.table(table(sample(beh_field$pollinatorClass, replace = TRUE)))
}

reps <- t(replicate(10000, resTab()))
colMeans(as.matrix(reps))
# get 95% CI's
propTab <- cbind(propTab, t(apply(reps, MARGIN = 2, quantile, c(0.025, 0.975))))
propTab
names(propTab)[3:4] <- c("lower", "upper")

for(ii in 1:6){
    print(prop.test(table(beh_field$pollinator)[ii], nrow(beh_field)))
}

nrow(beh_field)
ggplot(propTab, aes(x = Insect, y = Proportion)) + 
     #geom_bar(stat = 'identity') + 
     geom_point(size = 1, color= 'grey40') + 
     geom_pointrange(aes(ymax = upper, ymin=lower), size = 0.2) + 
     labs (x = "Insect visitor", y = "Proportion of visits") +
     ylim(c(0,1)) + 
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
           axis.ticks.x=element_blank()) + 
     scale_x_discrete(labels=c("Bombus spp." = expression(
          paste(italic("Bombus spp."))), 
          "Xylocopa virginica" = expression(
               paste(italic("Xylocopa virginica"))),
          "Apis mellifera" = expression(
               paste(italic("Apis mellifera")))))
ggsave(filename = "KalmiaFigures/VisitorProps2.pdf", width = 5, height = 4)

################################################################################
# look only at pollinators that triggered anthers
beh_trig <- beh_field[!(is.na(beh_field$timeTrig1)), ]
beh_trig <- droplevels(beh_trig)

propTab2 <- data.frame(prop.table(table(beh_trig$pollinatorClass)))
names(propTab2) = c("Insect", "Proportion")
propTab2

# resample proportion table to get error bars
resTab2 <- function(o){
     prop.table(table(sample(beh_trig$pollinatorClass, replace = TRUE)))
}

reps <- t(replicate(10000, resTab2()))
colMeans(as.matrix(reps))
# get 95% CI's
propTab2 <- cbind(propTab2, t(apply(reps, MARGIN = 2, quantile, c(0.025, 0.975))))
propTab2
names(propTab2)[3:4] <- c("lower", "upper")

ggplot(propTab2, aes(x = Insect, y = Proportion)) + 
     geom_point(size = 1, color= 'grey40') + 
     geom_pointrange(aes(ymax = upper, ymin=lower), size = 0.2) + 
     labs (x = "Insect visitor", y = "Proportion of visitors\nthat triggered anthers") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
           axis.ticks.x=element_blank()) + 
     scale_x_discrete(labels=c("Bombus spp." = expression(
          paste(italic("Bombus spp."))), 
          "Xylocopa virginica" = expression(
               paste(italic("Xylocopa virginica"))),
          "Apis mellifera" = expression(
               paste(italic("Apis mellifera")))))
ggsave(filename = "KalmiaFigures/VisitorTrig.pdf", width = 4, height = 4)
nrow(beh_trig)


######################################################################
## Describe what bumblebee behaviors resulted in anther triggers
######################################################################

# make a long dataset
newDF <- data.frame()
for(ii in beh$Video){
     tmp  = beh[beh$Video == ii , ]
     tmp_long <- reshape(tmp, varying = c("timeTrig1"    , "trigBehavior1" , "trigLeg1"     ,
                                          "timeTrig2",      "trigBehavior2",  "trigLeg2"   ,
                                          "timeTrig3"   ,  "trigBehavior3" , "trigLeg3"   ,
                                          "timeTrig4" ,     "trigBehavior4" , "trigLeg4"   ,
                                          "timeTrig5" ,    "trigBehavior5" , "trigLeg5"), 
             direction = 'long', 
             idvar = "Video", sep = "")
     
     tmp_long <- tmp_long[!(is.na(tmp_long$timeTrig)), ]
     newDF <- rbind(tmp_long, newDF)
}

# make new factors
newDF$Beh2 <- mapvalues(newDF$trigBehavior, 
                        from = c("inserting proboscis", 
                                "flapping wings", 
                                "walking" ,
                                "drinking", 
                                "collecting pollen" , 
                                'grooming'), 
                        to = c("proboscis extended", 
                               "wings flapping", 
                               "wings not flapping, proboscis retracted", 
                               "proboscis extended", 
                               "wings not flapping, proboscis retracted", 
                               "wings not flapping, proboscis retracted"))

# use only field-observations of bumblebees
newDF <- newDF[newDF$vid.type == 'field' & newDF$pollinator == 'bumblebee', ]

# redefine legs to just front, mid, rear, and unknown
newDF <- within(newDF, 
                trL2 <- mapvalues(trigLeg, from = c('right front', 'left front', 
                                                    'left mid', 'right mid'), 
                                  to = c('front leg', 'front leg', 'mid leg', 'mid leg')))


# set order of levels
newDF <- within(newDF, 
                trigBehavior <- factor(trigBehavior, 
                                   levels=names(sort(table(trigBehavior), 
                                                     decreasing=TRUE))))

newDF <- within(newDF, 
                Beh2 <- factor(Beh2, 
                                       levels=names(sort(table(Beh2), 
                                                         decreasing=TRUE))))
colnames(newDF)
ggplot(newDF, aes(x = trigBehavior)) + 
     geom_bar()

bb <- ggplot(newDF, aes(x = Beh2)) + 
     geom_bar() + 
     labs(x = "Behavior during pollen release") +  
     theme(legend.position = 'none') + 
     scale_fill_viridis(option = 'viridis', discrete = TRUE)
bb

ggplot(newDF, aes(x = 'Behavior that releases anther', fill = trigBehavior)) + 
     geom_bar() + 
     scale_fill_viridis(option = 'viridis', discrete = TRUE)



### bootstrap to get CI's (make a plot like bb)
vidnames = unique(newDF$Video)
behResamp <- function(o){
     # resample videos
     tmp = newDF[newDF$Video %in% sample(vidnames, replace = TRUE), c("Video", "Beh2") ]
     
     # resample behaviors within each video (if there's more than one)
     for(ii in names(table(tmp$Video))[table(tmp$Video) > 1]){
          tmp2 = tmp[tmp$Video == ii , ]
          tmp[tmp$Video == ii , ] =  tmp2[sample(1:nrow(tmp2), replace = TRUE), ]
     }
     
     pt = prop.table(table(tmp$Beh2))
     return(pt)
}


# should take about 30 seconds
system.time({
     beh_rep <- t(replicate(n = 10000, behResamp()))
})

colMeans(as.matrix(beh_rep))
# get 95% CI's
propTab3 <- data.frame(cbind(prop.table(table(newDF$Beh2)), 
                  t(apply(beh_rep, MARGIN = 2, 
                          quantile, c(0.025, 0.975)))))
colnames(propTab3) <- c('mean', "lower", "upper")
propTab3$Behavior <- rownames(propTab3)
propTab3


ggplot(propTab3, aes(x = Behavior, y = mean)) + 
     #geom_point(size = 1, color= 'grey40') + 
     geom_pointrange(aes(ymax = upper, ymin=lower), size = 0.2) + 
     labs (x = "Behavior", y = "Proportion of bumblebee behaviors\nthat triggered anthers") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
           axis.ticks.x=element_blank()) + 
     scale_x_discrete(labels=c("wings not flapping, proboscis retracted" = expression(
               paste("wings not flapping,\nproboscis retracted"))))
ggsave(filename = "KalmiaFigures/BehaviorTrig.pdf", width = 4, height = 4)
length(unique(newDF$Video))
nrow(newDF)



######################################################################
# Describe what leg triggered the anther (in bumblebees)
######################################################################


unique(newDF$trigLeg)

newDF <- within(newDF, 
                trigLeg <- factor(trigLeg, 
                               levels=names(sort(table(trigLeg), 
                                                 decreasing=TRUE))))

ggplot(newDF, aes(x = trigLeg, fill = vid.type)) + 
     geom_bar() + 
     labs(x = "What triggered the anther") +
     scale_fill_viridis(option = 'viridis', discrete = TRUE)



ggplot(newDF, aes(x = trL2, fill = vid.type)) + 
     geom_bar() + 
     labs(x = "What triggered the anther") + 
     scale_fill_viridis(option = 'viridis', discrete = TRUE)

# calculate CI's for proportions
### bootstrap to get CI's (make a plot like bb)
vidnames = unique(newDF$Video)
legResamp <- function(o){
     # resample videos
     tmp = newDF[newDF$Video %in% sample(vidnames, replace = TRUE), c("Video", "trL2") ]
     
     # resample behaviors within each video (if there's more than one)
     for(ii in names(table(tmp$Video))[table(tmp$Video) > 1]){
          tmp2 = tmp[tmp$Video == ii , ]
          tmp[tmp$Video == ii , ] =  tmp2[sample(1:nrow(tmp2), replace = TRUE), ]
     }
     
     pt = prop.table(table(tmp$trL2))
     return(pt)
}


# should take about 30 seconds
system.time({
     leg_rep <- t(replicate(n = 10000, legResamp()))
})

colMeans(as.matrix(leg_rep))
# get 95% CI's
propTab4 <- data.frame(cbind(prop.table(table(newDF$trL2)), 
                             t(apply(leg_rep, MARGIN = 2, 
                                     quantile, c(0.025, 0.975)))))
colnames(propTab4) <- c('mean', "lower", "upper")
propTab4$leg <- rownames(propTab4)
propTab4

propTab4$leg <- factor(propTab4$leg, levels = c('unknown', 'front leg', 'mid leg'))

### HERE

ggplot(propTab4, aes(x = leg, y = mean)) + 
     #geom_point(size = 1, color= 'grey40') + 
     geom_pointrange(aes(ymax = upper, ymin=lower), size = 0.2) + 
     labs (x = "Cause of anther release", 
           y = "Proportion of bumblebee causes\nthat triggered anthers") +
     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), 
           axis.ticks.x=element_blank()) + 
     ylim(c(0, 0.8))
ggsave(filename = "KalmiaFigures/LegTrig.pdf", width = 3, height = 4)
length(unique(newDF$Video))
nrow(newDF)
table(newDF$pollinator)

