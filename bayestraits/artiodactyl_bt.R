#### Artiodactyl BT ####

## packages
library(tidyverse)
library(btw)

## bt path
options(bt_path = ".",
        bt_bin = "BayesTraitsV3")

## data
dat <- read.table("bayestraits/Artiodactyl.txt")
tt <- ape::read.nexus("bayestraits/Artiodactyl.trees") ## 500 trees

## MCMC with exp 10 prior, default length and burnin
command_vec<-c("1","2","PriorAll exp 10")
mcmc_run<-bayestraits(data = dat, tree = tt, commands = command_vec)

## save or read in run
#saveRDS(mcmc_run, file = "bayestraits/artio_bt_mcmc.RDS")
mcmc_run<-readRDS("bayestraits/artio_bt_mcmc.RDS")

## ML for each tree
command_vec_ml<-c("1","1")
ml_run<-bayestraits(data = dat, tree = tt, commands = command_vec_ml)

## save or read in run
#saveRDS(ml_run, file = "bayestraits/artio_bt_ml.RDS")
ml_run<-readRDS("bayestraits/artio_bt_ml.RDS")

## get results
mcmc_results<-mcmc_run$Log$results
ml_results<-ml_run$Log$results

## plot rates?
mcmc_results %>% select(starts_with("q")) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x=name, y=value))+
  geom_boxplot()+
  theme_bw()

ml_results %>% select(starts_with("q")) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x=name, y=value))+
  geom_boxplot()+
  theme_bw()

## or to match with figure in manual, V3 page 21, looks pretty similar
mcmc_results %>% select(starts_with("q")) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x=value))+
  geom_histogram(binwidth = 1)+
  theme_bw()+
  facet_wrap(~name)

## parameter restriction example
mcmc_ind<-bayestraits(data = dat, tree = tt,
                      commands = c("1","2","PriorAll exp 10","Stones 100 1000"))
mcmc_dep<-bayestraits(data = dat, tree = tt,
                      commands = c("1","2","PriorAll exp 10","Stones 100 1000",
                                   "Restrict qDG qGD"))

## save/read results
#saveRDS(mcmc_ind, file = "bayestraits/artio_bt_mcmc_ind.RDS")
#saveRDS(mcmc_dep, file = "bayestraits/artio_bt_mcmc_dep.RDS")
mcmc_ind<-readRDS("bayestraits/artio_bt_mcmc_ind.RDS")
mcmc_dep<-readRDS("bayestraits/artio_bt_mcmc_dep.RDS")

## get BF
loglik_ind<-mcmc_ind$Stones$logMarLH # -8.775 in manual, -8.715 here
loglik_dep<-mcmc_dep$Stones$logMarLH # -8.296 in manual, -8.270 here
bf<-2*(loglik_ind-loglik_dep) # -0.957 in manual, -0.889 here
