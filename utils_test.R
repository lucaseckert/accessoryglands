library(lattice)
library(coda)
library(MCMCpack)
library(ramcmc)

## example from MCMCpack::MCMCmetrop1R
logitfun <- function(beta, y, X) {
  return(sum(dbinom(y, prob = plogis(X %*% beta), size = 1, log = TRUE)))
}
n <- 1000
set.seed(101)
X0 <- cbind(1, x1 = rnorm(n), x2 = rnorm(n))
beta0 <- c(0.5, -1, 1)
theta0 <- c(0, 0, 0)
y0 <- rbinom(n, size = 1, prob = plogis(X0 %*% beta0))
library(MCMCpack)
post.samp <- MCMCmetrop1R(logitfun, theta.init = theta0,
                          X = X0, y = y0,
                          thin = 1, mcmc = 40000, burnin = 500,
                          tune = c(1.5, 1.5, 1.5),
                          verbose = 500, logfun = TRUE)

xyplot(post.samp)

##
lfun2 <- function(beta, y, X) {
  return(sum(dbinom(y, prob = plogis(X %*% beta), size = 1, log = TRUE)))
}

opt1 <- optim(par = theta0, fn = logitfun, control = list(fnscale = -1),
              hessian = TRUE, X = X0, y = y0)
S <- t(chol(solve(-opt1$hessian)))
m <- metropolis(logitfun, theta0 = opt1$par, S = S,
                n_iter = 40000, n_burnin = 2000,
                p_args = list(X = X0, y = y0),
                adapt = TRUE)

xyplot(mcmc(m$theta))
m$S

sfun <- function() rnorm(3)
debug(metrop_mult)
m2 <- metrop_mult(logitfun, sfun, S = S,
                  nchains = 3, nclust = 3,
                  p_args = list(X = X0, y = y0),
                  n_iter=10000, n_burnin=500, adapt=TRUE)

xyplot(m2)
