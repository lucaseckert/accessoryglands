L <- load("cache/ag_fit.rda")
source("R/utils.R")
library(corHMM)

data <- ag_compdata$data
data0 <- ag_compdata$data0
CorData <- corHMM:::corProcessData(data, collapse = TRUE)
## pull machinery from corProcessData ?
tt <- traitNames(data0)   ## exclude first column (species names)
CorData$PossibleTraits

image(ag_model0, data0)  ## add numbers ?
m1 <- ag_model0$index.mat
dimnames(m1) <- list(tt, tt)

parnums <- c(na.omit(unique(c(m1))))
for (i in parnums) {
    w <- which(m1==i, arr.ind=TRUE, useNames=TRUE)
    w_to <- rownames(w)
    w_from <- colnames(m1)[w[, "col"]]
    ## figure out focal trait (split & compare)
    ss <- strsplit(rownames(w),"_")
    ssdf <- do.call(rbind,ss)
    labs <- apply(ssdf,2,
                  function(x) {
                      ## FIXME: > binary traits?
                      const <- (length(u <- unique(x)) ==1)
                      if (const) return(u)
                      return("")
                      ## return(gsub("[0-9]+$","",x))
                  })
    stopifnot(nrow(labs)==1) ## should be unique
    cat(i, paste(labs[nzchar(labs)], collapse="_"),"\n")
    ## still need to work out which trait is 'focal' changing (if more than one)
    ## "focal"
    ##
}
