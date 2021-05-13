#loading package
library(fishtree)

#loading data
data<-read.csv("binaryTraitData.csv",header = T)

#accessing tree block
treeBlock<-fishtree_complete_phylogeny(data$species)

#accessing single tree
treeSingle<-fishtree_phylogeny(data$species)

#saving
saveRDS(treeBlock, file = "treeBlock.rds")
saveRDS(treeSingle, file = "treeSingle.rds")
