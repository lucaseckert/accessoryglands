
load("cache/MK_3state_simple.rda")
library("MCMCpack")
logpostfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  loglik <- -1 * nllfun(p)
  log.prior <- sum(dnorm(p, mean = prior.mean, sd = prior.sd, log = TRUE))
  return(loglik + log.prior) ## product of likelihood and prior -> sum of LL and log-prior
}
m1 <- MCMCmetrop1R(logpostfun, p, verbose = 1000, thin=100, mcmc = 500000)
saveRDS(m1, file = "cache/MK_3state_constr1_mcmc.rds")
