dat <- read.table("bayestraits/Artiodactyl.txt")
tt <- ape::read.nexus("bayestraits/Artiodactyl.trees") ## 500 trees

system("BayesTraitsV4.1.2-Linux/BayesTraitsV4 bayestraits/Artiodactyl.trees bayestraits/Artiodactyl.txt < bayestraits/ArtiodactylMLIn.txt")
## runs ML on every tree ...


## scale trees

## set priors

##

## ML fits
corHMM(phy = ag_compdata_tb$phy,
       data = ag_compdata_tb$data,
       rate.cat = 1,
       rate.mat = ag_statemat_pcsc,
       root.p = root.p,
       lower = 0.1,                             ## 0.1 transitions per tree
       upper = 100 * ape::Ntip(ag_compdata$phy) ## 100 transitions per species
       )

## augment_model()?

## corhmm_mcmc()
