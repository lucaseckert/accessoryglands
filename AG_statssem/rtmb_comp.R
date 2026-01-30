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
