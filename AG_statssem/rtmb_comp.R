library(targets)
library(tarchetypes)
tar_load("ag_compdata")
tar_load("ag_statemat_pcsc")
tar_load("root.p")
tar_load("ag_corhmm_bounds")
source("R/utils.R")
library(corHMM)
statemat <- ag_statemat_pcsc
t1 <- system.time(
  a1 <- augment_model(
  corHMM(phy = ag_compdata$phy,
         data = ag_compdata$data,
         rate.cat = 1,
         rate.mat = statemat,
         root.p = root.p,
         lower.bound = ag_corhmm_bounds[["lower"]],
         upper.bound = ag_corhmm_bounds[["upper"]]
         )
))

t2 <- system.time(
  a2 <- update(a1, use_RTMB = TRUE)
)

a3 <- update(a1, use_RTMB = TRUE, return.devfun = TRUE)
all.equal(a1, a2, tolerance = 2e-3, check.attributes = FALSE)


## to make sure update() isn't doing some clever short-circuit?
## (e.g. starting from previous fitted parameter values?)
t4 <- system.time(
  a4 <- augment_model(
  corHMM(phy = ag_compdata$phy,
         data = ag_compdata$data,
         rate.cat = 1,
         rate.mat = statemat,
         root.p = root.p,
         lower.bound = ag_corhmm_bounds[["lower"]],
         upper.bound = ag_corhmm_bounds[["upper"]],
         use_RTMB = TRUE
         )
))

a4$opt.time
##   user  system elapsed 
##  3.943   0.001   3.943 
a1$opt.time
##    user  system elapsed 
##  56.214   0.009  56.220 

## 'only' a 15-fold speedup ...
## other components are taking some time ... what's slow?
## is it corHMM or augment.model?
## 31 vs 83 seconds ...
## computing Hessian -- could be sped up with RTMB as well ...

## make_nllfun
