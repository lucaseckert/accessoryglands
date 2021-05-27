do_MCMCpack <- FALSE
source("utils.R")

L <- load("cache/MK_3state_simple.rda")

logpostfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  loglik <- -1 * nllfun(p)
  log.prior <- sum(dnorm(p, mean = prior.mean, sd = prior.sd, log = TRUE))
  return(loglik + log.prior) ## product of likelihood and prior -> sum of LL and log-prior
}
assign("nllfun", nllfun, environment(logpostfun))

if (do_MCMCpack) {

  ## old: run single long chain via MCMC
  library("MCMCpack")
  m1 <- MCMCmetrop1R(logpostfun, p, verbose = 1000, thin=100, mcmc = 500000)
  saveRDS(m1, file = "cache/MK_3state_constr1_mcmc.rds")

} else {
  ## home-grown multi-chain adaptive MCMC
  library(bbmle)
  library(Matrix)
  ## need to "pos-defify" the cov matrix ...
  V <- as.matrix(nearPD(vcov(m0))$mat)
  max(eigen(V)$values)
  min(eigen(V)$values)
  S <- t(chol(V))
  make_sfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
    function() {
      prior.mean <- (lb + ub) / 2
      prior.sd <- (ub - lb) / (2 * range)
      rnorm(length(p), mean = prior.mean, sd = prior.sd)
    }
  }

  ## S2 <- diag(diag(S))
  S3 <- diag(rep(0.1, 12))  ## start with a sensible diagonal,
                            ## let adapt=TRUE adjust it
  set.seed(101)
  m1 <- metrop_mult(logpostfun,
              start = make_sfun(coef(m0), range=6),
              S = S3,
              nchains = 8, nclust = 8,
              n_burnin = 4000, n_iter = 44000, thin = 10,
              adapt = TRUE,
              clust_extra = lme4:::namedList(nllfun,
                                             MK_3state_simple))
  saveRDS(m1, file = "cache/MK_3state_constr1_mcmc_adapt.rds")

}
