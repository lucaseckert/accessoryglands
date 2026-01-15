dat <- read.table("bayestraits/Artiodactyl.txt")
tt <- ape::read.nexus("bayestraits/Artiodactyl.trees") ## 500 trees
remotes::install_github("bbolker/corHMM")
source("R/mcmc.R")
source("R/utils.R")
library(corHMM)

system("BayesTraitsV4.1.2-Linux/BayesTraitsV4 bayestraits/Artiodactyl.trees bayestraits/Artiodactyl.txt < bayestraits/ArtiodactylMLIn.txt")
## runs ML on every tree ...
## could read MLOut.txt

## scale trees
## "if no parameters are supplied the branches are scaled to have a mean of 0.1"
scaletree <- function(x, scale = 0.1) {
    x$edge.length <- scale*x$edge.length/mean(x$edge.length)
    return(x)
}
t0 <- scaletree(tt[[1]])
stopifnot(all.equal(mean(t0$edge.length), 0.1))

## example is NOT scaled
## tt <- lapply(tt, scaletree)

## set priors (?)

## how does BT handle rooting by default, and how do we mimic that?

## do an ML fit for a single tree
fitfun <- function(tree) {
    cat(".\n")  ## poor man's progress bar
    cc <- capture.output(res <- corHMM(phy = tree,
           data = dat,
           rate.cat = 1,
           rate.mat = matrix(1:4, 2,2),
           root.p = c(0.5, 0.5),
           ## "yang" doesn't work-- perhaps problem with only 2-state model?
           ## root.p = root.p,
           lower = 0.1,                             ## 0.1 transitions per tree
           upper = 100 * ape::Ntip(t) ## 100 transitions per species
           ))
    return(res)
}

fitfun(tt[[1]])

## ML fits
fn <- "artiodactyl_our_MLfits.rds"
if (!file.exists(fn)) {
    allfits <- lapply(tt, fitfun)
    saveRDS(allfits, fn)
}
## compare estimates with BT estimates ...

## augment_model()?

targets::tar_load("ag_model_pcsc")
ag_model_pcsc$args.list$p

## corhmm_mcmc()
pfun <- make_nllfun(fitfun(tt[[1]]), treeblock = tt)
corhmm_mcmc(fitfun(tt[[1]]),
            ## archaeology: not clear why we need nllfun in p_args?
            p_args = list(nllfun = make_nllfun(tt[[1]], treeblock = tt)),
            ## lb = log(1),
            ## ub = log(10 * ape::Ntip(tt[[1]]))),
            n_cores = 1,
            n_chains = 1,
            n_burnin = 4000,
            n_iter = 80000,
            n_thin = 10,
            seed = 101)
