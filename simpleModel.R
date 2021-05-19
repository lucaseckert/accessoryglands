#loading packages
if (packageVersion("corHMM")<"2.6.1") {
  stop("need hacked version of corHMM. Try 'remotes::install_github(\"bbolker/corHMM\", build_vignettes=TRUE)'")
}
## if installation fails when building vignettes drop build_vignettes=TRUE (default is FALSE)
library(corHMM)
library(bbmle)     ## for likelihood profiles/CIs
library(numDeriv)  ## for Hessian/Wald CIs
library(MCMCpack)  ## for Bayes (Metropolis-Hastings) run
library(coda)
library(lattice)
library(ggplot2); theme_set(theme_bw())
library(colorspace)
library(GGally)
library(ggthemes)
library(fishtree)
library(caper)

#loading data and tree and trimming
allData<-read.csv("binaryTraitData.csv",header = T)
fullPhy<-fishtree_phylogeny(allData$species)
traitData<-data.frame(species=allData$species, ag=allData$ag, care=allData$care, mating=allData$mating, names=allData$species)
trimmedData<-comparative.data(fullPhy,traitData, names.col = names)
phy<-trimmedData$phy
data<-trimmedData$data

#default model
(MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1))

#getting state matrix and legend
getStateMat4Dat(data)
LegendAndRateMat<-getStateMat4Dat(data)
RateMat<-LegendAndRateMat$rate.mat
RateMat

#making the gain and loss of care and sperm comp exogenous, down to 12 pars
pars2equal<-list(c(7,10,20,23), c(4,11,17,24), c(2,5,15,18), c(1,8,14,21))
StateMatA_constrained<-equateStateMatPars(RateMat, pars2equal)
StateMatA_constrained

#making the simplified model
MK_3state_simple <- corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = StateMatA_constrained)
MK_3state_simple

#I'm not sure how to get your hacked version of the package so this is as far as I can get 
