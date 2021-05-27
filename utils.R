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
#' @param trace_level verbosity (0=none; 1=log-posterior + accept rate; 2 includes parameter values)
#' @param trace_steps frequency of reports
#' @param seed random number seed
#' @param p_args list of additional arguments to pass to postfun
## FIXME: add trace/verbose; auto-optim? return as mcmc object + attributes? multi-chain version?
metropolis <- function(postfun,
                       theta0,
                       S,
                       n_iter=1000,
                       n_burnin=100,
                       thin = 1,
                       adapt = FALSE,
                       p_args = list(),
                       trace_level = 1,
                       trace_steps = 100,
                       seed = NULL) {
  if (adapt && !require("ramcmc")) {
    stop("need ramcmc package for adaptive MCMC")
  }

  ## FIXME: seed for parallel versions?
  p <- length(theta0)
  theta <- matrix(NA, nrow = ceiling((n_iter - n_burnin) / thin), ncol = p)
  accept <- numeric(n_iter)
  bad_steps <- 0

  pwrap <- function(x) {
    do.call("postfun",c(list(x), p_args))
  }

  ## initialize
  S0 <- S   ## store initial cholesky factor (for comparison later if we need)
  posterior <- pwrap(theta0)
  theta_current <- theta0
  j <- 0  ## storage counter

  for (i in 2:n_iter) {
    u <- rnorm(p)
    theta_prop <- theta_current + S %*% u
    posterior_prop <- pwrap(theta_prop)
    acceptance_prob <- min(1, exp(posterior_prop - posterior))
    if (!is.finite(acceptance_prob)) {
      bad_steps <- bad_steps + 1
      acceptance_prob <- 0
    }
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
    if (trace_level > 0 && i %% trace_steps == 0) {
      cat(sprintf("** %d: log_post = %1.2f accept_rate=%1.2f\n",
                  i, posterior, mean(accept[1:i])))
      if (trace_level > 1) {
        cat(theta_current,"\n")
      }
    }
  }
  list(theta = theta, S = S,
       acceptance_rate = sum(accept[(n_burnin + 1):n_iter]) / (n_iter - n_burnin))
}


#' @param postfun log-posterior
#' @param start randomized starting function, or list of starting values
#' @param S initial Cholesky factor for candidate distribution
#' @param nchains number of chains
#' @param nclust number of cores
#' @param clust_extras additional objects needed in clust environments
#' @param ... additional arguments to \code{metropolis}
metrop_mult <- function(postfun, start, S, nchains,
                        nclust = min(nchains,getOption("mc.cores", 1)),
                        clust_extras = list(),
                        ...) {

  ## FIXME: better/less clunky way to parallelize?
  ## (1) handling environments/objects; (2) progress bar?
  ## (3) random-number-seed handling
  if (is.function(start)) {
    start <- replicate(nchains, start(), simplify=FALSE)
  }
  if (length(start) != nchains) stop("need as many start values as chains")
  require("parallel")
  mfun <- function(s) {
    metropolis(postfun, theta0 = s, S = S, ...)
  }
  if (nclust == 1) {
    res <- lapply(start, mfun)
  } else {
    clust_stuff <- c(list(...), clust_extras)
    suppressMessages(
        attach(clust_stuff) ## yes, we really want/need to do this
    )
    on.exit(detach(clust_stuff))
    cl <- makeCluster(nclust)
    clusterExport(cl, varlist = c("mfun", "metropolis", "postfun", "S",
                                  names(clust_stuff)),
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
