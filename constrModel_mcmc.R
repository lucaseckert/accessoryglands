do_MCMCpack <- FALSE

L <- load("cache/MK_3state_simple.rda")

logpostfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  loglik <- -1 * nllfun(p)
  log.prior <- sum(dnorm(p, mean = prior.mean, sd = prior.sd, log = TRUE))
  return(loglik + log.prior) ## product of likelihood and prior -> sum of LL and log-prior
}

if (do_MCMCpack) {
  library("MCMCpack")
  source("utils.R")
  m1 <- MCMCmetrop1R(logpostfun, p, verbose = 1000, thin=100, mcmc = 500000)
  saveRDS(m1, file = "cache/MK_3state_constr1_mcmc.rds")
} else {
  ## home-grown multi-chain adaptive MCMC
  library(bbmle)
  library(Matrix)
  ## need to "pos-defify" the cov matrix ...
  V <- as.matrix(nearPD(vcov(m0))$mat)
  S <- t(chol(V))
  make_sfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
    function() {
      prior.mean <- (lb + ub) / 2
      prior.sd <- (ub - lb) / (2 * range)
      rnorm(length(p), mean = prior.mean, sd = prior.sd)
    }
  }

  metrop_mult(logpostfun,
              start = make_sfun(coef(m0)),
              S = S,
              nchains = 4, nclust = 1,
              n_burnin = 1000, n_iter = 2000,
              adapt=TRUE)

}


