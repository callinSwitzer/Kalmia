file.choose()

fls = dir("/Users/callinswitzer/Desktop/FastecSummer2016KalmiaVideos")
fls = fls[grep('.avi', fls)]

fls = fls[order(fls, decreasing = FALSE)]

fls = sapply(fls, FUN = function(x) strsplit(x, split = "\\.")[[1]][1])


write.csv(fls, file = "~/Desktop/HSVID.csv", row.names = FALSE)
