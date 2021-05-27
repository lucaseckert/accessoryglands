#' utility function for hexbin panels for MCMC pairs plot
#' @param data ...
#' @param mapping ...
#' @param lb lower bound
#' @param ub upper bound
#' @param range distance from center to bound, in SD
#' @param show_prior include representation of prior (red lines)
my_mcmc <- function(data, mapping, lb=log(1e-9), ub=log(1e2), range=2,
                    show_prior=FALSE, geom=c("hexbin","density"),
                    density_alpha=0.5,...) {
  geom <- match.arg(geom)
  gg1 <- ggplot(data = data, mapping=mapping)
  gg1 <- switch(geom,
                hexbin = gg1 + geom_hex()  + scale_fill_viridis_c(),
                density = gg1 + geom_density_2d_filled(alpha=density_alpha) + scale_fill_grey(start=0.9, end=0.02)
                )
    ## geom_point(..., alpha = 0.2)
  if (show_prior) {
    prior.mean <- (lb+ub)/2
    prior.sd <- (ub-lb)/(2*range)
    gg1 <- (gg1
      + geom_vline(xintercept=c(prior.mean+c(-2,0,2)*prior.sd),
                   lty=c(2,1,2),
                   colour="red",
                   alpha=0.3)
      + geom_hline(yintercept=c(prior.mean+c(-2,0,2)*prior.sd),
                   lty=c(2,1,2),
                   colour="red",
                   alpha=0.3)
    )
  }
  return(gg1)
}


#' adaptive Metropolis, modified from ramcmc package example
#' @param postfun log-posterior  probability function
#' @param theta0 starting value
#' @param S initial Cholesky factor for MVN candidate distribution
#' @param n_iter number of iterations
#' @param n_burnin length of burnin/adapt phase
#' @param thin thinning frequency
#' @param adapt (logical) adapt?
#' @param ... additional arguments to pass to postfun
## FIXME: add trace/verbose; auto-optim? return as mcmc object + attributes? multi-chain version?
metropolis <- function(postfun, theta0, S, n_iter=1000,
                       n_burnin=100, thin = 1, adapt = FALSE, p_args=list()) {
  if (adapt && !require("ramcmc")) {
    stop("need ramcmc package for adaptive MCMC")
  }
  ## FIXME: test/check; implement thinning; wrap in parallel version?
  p <- length(theta0)
  theta <- matrix(NA, nrow = ceiling((n_iter - n_burnin) / thin), ncol = p)
  accept <- numeric(n_iter)

  pwrap <- function(x) {
    do.call("postfun",c(list(x), p_args))
  }

  ## initialize
  posterior <- pwrap(theta0)
  theta_current <- theta0
  j <- 0  ## storage counter

  for (i in 2:n_iter) {
    u <- rnorm(p)
    theta_prop <- theta_current + S %*% u
    posterior_prop <- pwrap(theta_prop)
    acceptance_prob <- min(1, exp(posterior_prop - posterior))
    if (runif(1) < acceptance_prob) {
      accept[i] <- 1
      theta_current <- theta_prop
      posterior <- posterior_prop
    }
    if (i > n_burnin && i %% thin == 0) {
      j <- j + 1
      theta[j, ] <- theta_current
    }
    if (adapt && i <= n_burnin) {
      S <- ramcmc::adapt_S(S, u, acceptance_prob, i - 1)
    }
  }
  list(theta = theta, S = S,
       acceptance_rate = sum(accept[(n_burnin + 1):n_iter]) / (n_iter - n_burnin))
}


metrop_mult <- function(postfun, start, S, nchains, nclust = getOption("mc.cores", 1), ...) {
  if (is.function(start)) {
    start <- replicate(nchains, start(), simplify=FALSE)
  }
  require("parallel")
  mfun <- function(s) {
    metropolis(postfun, theta0 = s, S = S, ...)
  }
  if (nclust == 1) {
    res <- lapply(start, mfun)
  } else {
    L <- list(...)
    attach(L) ## yes, we really want/need to do this
    on.exit(detach(L))
    cl <- makeCluster(nclust)
    clusterExport(cl, varlist = c("mfun", "metropolis", "postfun", "S",
                                  names(L)),
                  envir = environment(mfun))
    res <- parLapply(cl, start, mfun)
    stopCluster(cl)
  }
  mval <- coda::as.mcmc.list(lapply(res,
                                    function(x) coda::as.mcmc(x$theta)))
  attr(mval, "acceptance_rate") <- vapply(res, "[[", FUN.VALUE=numeric(1),
                                          "acceptance_rate")
  return(mval)
}
