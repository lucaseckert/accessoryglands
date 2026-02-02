library(targets)
library(tarchetypes)
tar_load("ag_compdata")
tar_load("ag_statemat_pcsc")
tar_load("ag_statemat_full")
tar_load("root.p")
tar_load("ag_corhmm_bounds")
source("R/utils.R")
library(corHMM)


t1 <- system.time(
  a1 <- augment_model(
  corHMM(phy = ag_compdata$phy,
         data = ag_compdata$data,
         rate.cat = 1,
         rate.mat = ag_statemat_pcsc,
         root.p = root.p,
         lower.bound = ag_corhmm_bounds[["lower"]],
         upper.bound = ag_corhmm_bounds[["upper"]]
         )
))

t2 <- system.time(
  a2 <- update(a1, use_RTMB = TRUE)
)

all.equal(a1, a2, tolerance = 2e-3, check.attributes = FALSE)
all.equal(logLik(a1), logLik(a2))
all.equal(logLik(a1), logLik(a2), scale = 1) ## absolute, not relative comparison, is more appropriate
logLik(a1) - logLik(a2)

a1$opt.time/a2$opt.time

## try 24-parameter model
a_full <- update(a1, rate.mat = ag_statemat_full,
                 use_RTMB = TRUE)

a_full_corHMM <- update(a1, rate.mat = ag_statemat_full)


a5 <- corHMM(phy = ag_compdata$phy,
       data = ag_compdata$data,
       rate.cat = 1,
       rate.mat = ag_statemat_pcsc,
       root.p = root.p,
       lower.bound = ag_corhmm_bounds[["lower"]],
       upper.bound = ag_corhmm_bounds[["upper"]],
       use_RTMB = TRUE
       )

a5r <-update(a5, return.devfun = TRUE)
str(a5r$devfun)
ff <- a5r$devfun()
ff$fn()
ff$gr()
library(tmbstan)
tt <- tmbstan(ff)


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
