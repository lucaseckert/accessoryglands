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


#' adaptive Metropolis from ramcmc package
#' @param postfun log-posterior  probability function
#' @param theta0 starting value
#' @param S initial Cholesky factor for MVN candidate distribution
#' @param n_iter number of iterations
#' @param n_burnin length of burnin/adapt phase
#' @param adapt (logical) adapt?
metropolis <- function(postfun, theta0, S, n_iter, n_burnin, adapt = FALSE) {

  ## FIXME: test/check; implement thinning; wrap in parallel version?
  p <- length(theta0)
  theta <- matrix(NA, n_iter, p)
  accept <- numeric(n_iter)
  posterior <- postfun(theta0)
  theta[1, ] <- theta0

  for (i in 2:n_iter){
    u <- rnorm(p)
    theta_prop <- theta[i - 1, ] + S %*% u
    posterior_prop <- postfun(theta_prob)
    acceptance_prob <- min(1, exp(posterior_prop - posterior))
    if (runif(1) < acceptance_prob) {
      accept[i] <- 1
      theta[i, ] <- theta_prop
      posterior <- posterior_prop
    } else{
      theta[i, ] <- theta[i - 1, ]
    }
    if(adapt && i <= n_burnin) {
      S <- ramcmc::adapt_S(S, u, acceptance_prob, i - 1)
    }
  }
  list(theta = theta[(n_burnin + 1):n_iter, ], S = S,
       acceptance_rate = sum(accept[(n_burnin + 1):n_iter]) / (n_iter - n_burnin))
}
