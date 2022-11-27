x <- read.csv("rjBayesTraits.csv")
rmat <- do.call("rbind",strsplit(x$Model.string, " "))
mean(rmat[,2] == rmat[,4])
