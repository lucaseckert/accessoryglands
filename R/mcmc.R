#' adaptive Metropolis, modified from ramcmc package example
#' @param postfun log-posterior  probability function
#' @param theta0 starting value
#' @param S initial Cholesky factor for MVN candidate distribution
#' @param n_iter number of iterations
#' @param n_burnin length of burnin/adapt phase
#' @param n_thin thinning frequency
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
                       n_thin = 1,
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
  theta <- matrix(NA, nrow = ceiling((n_iter - n_burnin) / n_thin), ncol = p)
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
    if (i > n_burnin && i %% n_thin == 0) {
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
#' @param n_chains number of chains
#' @param n_cores number of cores
#' @param clust_extras additional objects needed in clust environments
#' @param ... additional arguments to \code{metropolis}
metrop_mult <- function(postfun, start, S, n_chains,
                        n_cores = min(n_chains,getOption("mc.cores", 1)),
                        clust_extras = list(),
                        ...) {

  ## FIXME: better/less clunky way to parallelize?
  ## (1) handling environments/objects; (2) progress bar?
  ## (3) random-number-seed handling
  if (is.function(start)) {
    start <- replicate(n_chains, start(), simplify=FALSE)
  }
  if (length(start) != n_chains) stop("need as many start values as chains")
  require("parallel")
  mfun <- function(s) {
    metropolis(postfun, theta0 = s, S = S, ...)
  }
  if (n_cores == 1) {
    res <- lapply(start, mfun)
  } else {
    clust_stuff <- c(list(...), clust_extras)
    suppressMessages(
        attach(clust_stuff) ## yes, we really want/need to do this
    )
    on.exit(detach(clust_stuff))
    if (Sys.getenv("RSTUDIO") == "1" &&
        !nzchar(Sys.getenv("RSTUDIO_TERM"))
        &&  Sys.info()["sysname"] == "Darwin"
        && getRversion() >= "4.0.0") {
      cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
    } else {
      cl <- parallel::makeCluster(n_cores)
    }
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

##' @param model log-posterior function
##' @param n_cores number of parallel cores
##' @param n_chains
##' @param n_burnin
##' @param n_iter
##' @param n_thin
##' @param start_cand_range
##' @param seed
##' @param p_args additional arguments to pass to the log-posterior function
##' c1 <- corhmm_mcmc(ag_model0, p_args=list(nllfun=make_nllfun(ag_model0)), n_cores = 2, n_burnin=0, n_iter=100)
corhmm_mcmc <- function(model,
                        n_cores = NULL,
                        n_chains = 8,
                        n_burnin = 4000,
                        n_iter = 44000,
                        n_thin = 10,
                        start_cand_range = 6,
                        seed = NULL,
                        p_args = NULL,
                        treeblock = NULL) {

  if (is.null(n_cores)) {
    n_cores <- min(n_chains, getOption("mc.cores", 2))
  }

  require("bbmle")
  f <- make_nllfun(model, treeblock = treeblock)
  p <- model$args.list$p

  if (!is.null(seed)) set.seed(seed)

  if (FALSE) {
    ## attempt to set starting candidate distribution based on Hessian
     V <- solve(numDeriv::hessian(f, p))
     if (min(eigen(V)$values) < 1e-16) {
       V <- as.matrix(Matrix::nearPD(V)$mat)
     }
     S <- t(chol(V))
  }

  ## or: start with a sensible diagonal, let adapt=TRUE adjust it
  S <- diag(rep(0.1, length(p)))

  make_sfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
    function() {
      prior.mean <- (lb + ub) / 2
      prior.sd <- (ub - lb) / (2 * range)
      rnorm(length(p), mean = prior.mean, sd = prior.sd)
    }
  }

  m1 <- metrop_mult(
      corhmm_logpostfun,
      start = make_sfun(p, range=start_cand_range),
      S = S,
      n_chains = n_chains,
      n_cores = n_cores,
      n_burnin = n_burnin,
      n_iter = n_iter,
      n_thin = 10,
      adapt = TRUE,
      p_args = p_args,
      clust_extra = tibble::lst(model, corhmm_logpostfun, treeblock))

  m1[] <- lapply(m1, function(x) { dimnames(x) <- list(NULL, names(model$args.list$p)); return(x) })

  return(m1)
}
