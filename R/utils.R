##' utility function for hexbin panels for MCMC pairs plot
##' @param data ...
##' @param mapping ...
##' @param lb lower bound
##' @param ub upper bound
##' @param range distance from center to bound, in SD
##' @param show_prior include representation of prior (red lines)
my_mcmc <- function(data, mapping, lb=log(1e-9), ub=log(1e2), range=2,
                    show_prior=FALSE, geom=c("hexbin","density"),
                    density_alpha=0.5,
                    hpd_levels = c(0.5, 0.8, 0.9, 0.95),
                    ...) {
  require("emdbook") ## FIXME: not sure why? remove & see what breaks?
  geom <- match.arg(geom)
  gg1 <- ggplot(data = data, mapping=mapping)
  ## if (geom=="hexbin") browser()
  gg1 <- switch(geom,
                hexbin = gg1 + geom_hex()  + scale_fill_viridis_c(),
                ## inefficient because we are computing
                density = {
                  dd <- ggplot_build(gg1)$data[[1]]
                  levels <- get_hpd2d_levels(dd$x, dd$y, prob= hpd_levels)
                  (gg1
                    + geom_density_2d_filled(alpha=density_alpha)
                    + geom_density_2d(breaks = levels, colour="red")
                    + scale_fill_grey(start=0.9, end=0.02)
                  )
                }
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

##' Get 2D highest posterior density levels corresponding to probability regions
##' @param x x-coordinate of samples
##' @param y y-coordinate
##' @param probs vector of probability levels
##' @param ... arguments for MASS::kde2d
##' @examples
##' dd <- data.frame(x=rnorm(1000), y=rnorm(1000))
##' get_hpd2d_levels(dd$x,dd$y)
##' gg2 <- ggplot(dd, aes(x,y)) + geom_density_2d(breaks=get_hpd2d_levels(dd$x,dd$y), colour="red")
##' print(gg2)
get_hpd2d_levels <- function(x, y, prob=c(0.9,0.95), ...) {
  post1 <- MASS::kde2d(x, y)
  dx <- diff(post1$x[1:2])
  dy <- diff(post1$y[1:2])
  sz <- sort(post1$z)
  c1 <- cumsum(sz) * dx * dy
  ## remove duplicates
  ## dups <- duplicated(sz)
  ## sz <- sz[!dups]
  ## c1 <- c1[!dups]
  levels <- sapply(prob, function(x) {
    approx(c1, sz, xout = 1 - x, ties = mean)$y
  })
  return(levels)
}

##



#' @param dd trait data (no species names)
#' @param sep1
#' @param sep2
state_names <- function(dd, sep1 = "", sep2 = "_") {
  get_u <- function(x) {
    s <- sort(unique(x))
    s[s != "?"]
  }
  uu <- lapply(dd, get_u)
  uu2 <- mapply(paste, names(uu), uu, MoreArgs = list(sep = sep1),
                SIMPLIFY = FALSE)
  ## expand in slowest - to -fastest order
  uu3 <- rev(do.call("expand.grid", rev(uu2)))
  return(apply(uu3, 1, paste, collapse = sep2))
}

mk_idf <- function(index.mat) {
  df <- data.frame(which(!is.na(index.mat), arr.ind=TRUE),
                   index=c(na.omit(c(index.mat))))
  return(df)
}

## my.legend <- packGrob(
##     draw.colorkey(key=list(col = rainbow, at = do.breaks(range(z),100))), 
##     textGrob("My title", x = 0, y = 0.5, just = c("left", "centre")), 
##     height=unit(2, "lines"),side="top", dynamic=TRUE)

##' utility to plot transition matrix
##' @param x transition matrix
##' @param iso aspect ratio
plot.corhmm <- function(x,
                        aspect="iso",
                        log = TRUE,
                        include_nums = TRUE,
                        pnum_col = "red",
                        pnum_cex = 1.5,
                        ...) {
  require("Matrix")
  require("latticeExtra")
  dd <- x$data[, -1] ## no species names
  nm <- state_names(dd)
  M <- x$solution
  if (log) M <- log10(M)
  rlabs <- clabs <- nm
  ## YES, for the 1000th time "from ROW to COLUMN" (see corHMM vignette)
  p <- Matrix::image(Matrix(M),
                       scales=list(x=list(at=seq(nrow(M)),labels=rlabs,
                                          rot=90),
                                   y=list(at=seq(ncol(M)),labels=clabs)),
                       xlab="to",
                       ylab="from",
                       sub="",
                       aspect=aspect, ...)
  if (include_nums) {
      dd <- c(mk_idf(x$index.mat), pnum_col=pnum_col, pnum_cex=pnum_cex)
      ## n.b. col = from = x, row = to = y
      p <- p + latticeExtra::layer(lattice::panel.text(col, row, index, col=pnum_col, cex=pnum_cex), data=dd)
  }
  return(p)
}


iplot <- function(M, rlabs = rownames(M), clabs = colnames(M), aspect = "iso") {
    Matrix::image(Matrix(M),
                  scales=list(x=list(at=seq(nrow(M)),labels=rlabs,
                                     rot=90),
                              y=list(at=seq(ncol(M)),labels=clabs)),
                  xlab="to",
                  ylab="from",
                  sub="",
                  aspect=aspect)
}

##' wrap corhmm results to create a negative log-likelihood function
##' (faster than calling corHMM(p = exp(log_p))
##' @param corhmm_fit
##' @param treeblock sample a random phylogeny from the treeblock when computing?
make_nllfun <- function(corhmm_fit, treeblock = NULL) {
  require("bbmle")
  f <- function(log_p, treeblock = NULL) {
    a <- corhmm_fit$args.list
    a$p <- log_p
    if (!is.null(treeblock)) {
      a$phy <- treeblock[[sample(length(treeblock), size = 1)]]
    }
    return(do.call(corHMM:::dev.corhmm, a))
  }
  p <- corhmm_fit$args.list$p
  parnames(f) <- if (is.null(names(p))) {
                   paste0("p",seq_along(p))
                 } else names(p)
  return(f)
}

#' @param parameters (log-hazard rates)
#' @param lb lower bound(s) for baseline priors
#' @param ub upper bound(s)
#' @param range width of Gaussian (+/- SD between mean and lower/upper bounds)
#' @param gainloss_pairs
#' @param lb_gainloss
#' @param ub_gainloss
#' @param range_gainloss number of SDs from center to lower/upper bounds
#' @param nllfun \emph{negative} log-likelihood function
#' @param negative return negative log posterior?
corhmm_logpostfun <- function(p,
                              lb = log(1e-9),
                              ub = log(1e2),
                              range = 3,
                              gainloss_pairs = NULL,
                              lb_gainloss = log(1e-3),
                              ub_gainloss = log(1e3),
                              range_gainloss = 3,
                              nllfun,
                              negative = FALSE
                              ) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  loglik <- -1*nllfun(p)
  log.prior <- sum(dnorm(p, mean = prior.mean, sd = prior.sd, log = TRUE))
  ## product of likelihood and prior -> sum of LL and log-prior
  res <- loglik + log.prior
  ## calculate gain/loss priors
  if (!is.null(gainloss_pairs)) {
    gl.prior.mean <- (lb_gainloss + ub_gainloss) / 2
    gl.prior.sd <- (ub_gainloss - lb_gainloss) / (2 * range_gainloss)
    gl.values <- vapply(gainloss_pairs,
                        function(x) p[x[1]] - p[x[2]],
                        FUN.VALUE = numeric(1))
    gl.log.prior <- sum(dnorm(gl.values, mean = gl.prior.mean, sd = gl.prior.sd,
                              log = TRUE))
    res <- res + gl.log.prior
  }
  if (negative) res <- -1*res
  return(res)
}

#' @param model a \code{corHMM} model
#' @param data trait data associated with the model (matrix of traits only)
#' @examples
#' load("cache/ag_fit.rda")
#' par_names(ag_model0, ag_compdata$data[,-1])
par_names <- function(model) {
    data <- model$data[, -1]
    m1 <- model$index.mat
    tt <- state_names(data)
    dimnames(m1) <- list(to = tt, from = tt)
    parnums <- c(na.omit(unique(c(m1))))
    nm <- character(length(parnums))
    for (i in parnums) {
        w <- which(m1==i, arr.ind=TRUE, useNames=TRUE)
        w_to <- rownames(w)
        w_from <- colnames(m1)[w[, "from"]]
        ss_to <- strsplit(w_to, "_")
        ss_from <-strsplit(w_from, "_")
        focal <- ss_to[[1]][ss_to[[1]] != ss_from[[1]]]
        ssdf <- do.call(rbind,ss_to)
        labs <- apply(ssdf,2,
                      function(x) {
                          ## FIXME: > binary traits?
                          const <- (length(u <- unique(x)) ==1)
                          if (const) return(u)
                          return("")
                          ## return(gsub("[0-9]+$","",x))
                      })
        stopifnot(nrow(labs)==1) ## should be unique
        lab <- paste(labs[nzchar(labs)], collapse="_")
        lastnum <- as.numeric(substr(focal, nchar(focal), nchar(focal)))
        focal_lab <- paste(c("gain", "loss")[lastnum + 1],
                           substr(focal, 1, nchar(focal) - 1),
                           sep=".")
        lab <- gsub(focal, focal_lab, lab)
        nm[i] <- lab
    }
    return(nm)
}


#' give parameters and states more meaningful names
#' @param model corHMM model
#' @param add_hessian compute Hessian and add it as an attribute?
augment_model <- function(model, add_hessian = TRUE) {
    sn <- state_names(model$data[, -1])
    dimnames(model$solution) <- dimnames(model$index.mat) <- list(sn, sn)
    names(model$args.list$p) <- par_names(model)
    if (add_hessian) {
        ## update model with return.devfun?
        f <- make_nllfun(model)
        p <- model$args.list$p
        H <- numDeriv::hessian(f, p)
        ## if we had RTMB functions we could use numDeriv::jacobian(gr, p)
        attr(model, "hessian") <- H
    }
    return(model)
}

##' @inheritParams tidy.mle2
##' @param prof_args arguments to pass to profile.corhmm()
##' @param contrast_mat contrast matrix (right-multiply by parameters)
tidy.corhmm <- function(x,
                        conf.int = FALSE,
                        conf.method = "wald",
                        conf.level = 0.95,
                        p.value = FALSE,
                        profile = NULL,
                        exponentiate = FALSE,
                        prof_args = NULL,
                        contrast_mat = NULL,
                        ...) {
    f <- make_nllfun(x)
    ## save *original* p in case we use contrasts
    p0 <- p <- x$args.list$p

    if (!is.null(contrast_mat)) {
        p <- drop(p %*% contrast_mat)
    }
    ## FIXME: augment model here ? or assume already augmented?
    res <- dplyr::tibble(term = names(p),
                         estimate = p)
    if (conf.int) {
        if (conf.method == "wald") {
            ## FIXME: expensive to compute. Upstream? Cache?
            H <- attr(x, "hessian")
            if (is.null(H)) {
                H <- numDeriv::hessian(f, p0)
            }
          V <- solve(H)
          if (!is.null(contrast_mat)) {
              V <- t(contrast_mat) %*% V %*% contrast_mat
          }
          sds <- sqrt(diag(V))
          qq <- qnorm((1+conf.level)/2)
          res <- dplyr::mutate(res,
                    std.error = sds,
                    statistic = p/sds,
                    conf.low = estimate - qq*std.error,
                    conf.high = estimate + qq*std.error)
    } else if (conf.method == "profile") {
        if (is.null(profile)) {
            profile <- do.call(profile.corhmm, c(list(x), prof_args))
        }
        cfun <- function(x) {
            junk <- capture.output(r <- try(confint(x), silent=TRUE))
            if (inherits(r, "try-error")) return(data.frame(conf.low = NA, conf.high = NA))
            return(setNames(r, c("conf.low", "conf.high")))
        }
        cc <- suppressWarnings(purrr:::map_dfr(profile, cfun, .id="term"))
        res <- full_join(res, cc, by = "term")
    }
  }
  if (exponentiate) {
    res <- mutate(res, across(estimate, exp))
    res <- mutate(res, across(starts_with("conf"), exp))
    if ("std.error" %in% names(res)) {
      res <- mutate(res, across(std.error, ~ . * estimate))
    }
  }
  return(res)
}

glance.corhmm <- function(x, nobs = NULL, ...) {
    L <- logLik(x)
    res <-tibble(
        df = attr(L, "df"),
        logLik = c(L),
        AIC = 2*df - 2 * logLik)
    if (is.null(nobs)) {
        nobs <- tryCatch(nobs(x),
                         error = function(e) NULL)
    }
    if (!is.null(nobs)) {
        res <- mutate(res, BIC = -2*logLik + log(nobs)*df,
                   AICc = AIC + 2*(df^2 + df)/(nobs-df-1))
    }
    return(res)
}

glance.corhmm_contrast <- glance.mle2 <- glance.corhmm

profile.corhmm <- function(model,
                           mlefun = NULL,
                           which = NULL,
                           quietly = FALSE,
                           n_cores = getOption("mc.cores", 2), ...) {
  if (is.null(mlefun)) {
    mlefun <- make_nllfun(model)
  }
  p <- model$args.list$p
  ## FIXME:: use maxit=0 ?? (still need to compute Hessian,
  ##  unless prespecified)
  if (!quietly) cat("re-fitting with mle2\n")
  mlefun <- bbmle::mle2(mlefun, start = p)
  if (is.null(which)) which <- seq_along(p)
  pfun <- function(i) {
    profile(mlefun, which=i, ...)
  }
  if (n_cores==1) {
    res <- lapply(which, pfun)
  } else {
    ## FIXME: should build parallel capability into bbmle instead!
    if (!quietly) cat("executing on",n_cores,"cores\n")
    ## work around Rstudio/R bug
    ## https://stackoverflow.com/questions/62730783/error-in-makepsockclusternames-spec-cluster-setup-failed-3-of-3-work
    if (Sys.getenv("RSTUDIO") == "1" &&
        !nzchar(Sys.getenv("RSTUDIO_TERM"))
        &&  Sys.info()["sysname"] == "Darwin"
        && getRversion() >= "4.0.0") {
      cl <- parallel::makeCluster(n_cores, setup_strategy = "sequential")
    } else {
      cl <- parallel::makeCluster(n_cores)
    }
    parallel::clusterEvalQ(cl, library(bbmle))
    res <- parallel::parLapply(cl = cl, X = which, fun = pfun)
  }
  if (!quietly) cat("... done\n")
  ## could try to collapse into a single profile.mle2 object ...
  names(res) <- names(p)
  class(res) <- "profile.corhmm"
  return(res)
}

confint.profile.corhmm <- function(x) {
  ## indiv components are profile.mle2 objects
  require("bbmle")
  lapply(x, confint)
}

as.data.frame.profile.corhmm <- function(x) {
  ## indiv components are profile.mle2 objects
  require("bbmle")
  lapply(x, as.data.frame)
}

## copied from emdbook::lump.mcmc.list
as.mcmc.mcmc.list <- function(x) {
    x2 <- do.call("rbind", x)
    mcpars <- sapply(x, attr, "mcpar")
    class(x2) <- "mcmc"
    if (var(mcpars[1, ]) > 0 || var(mcpars[3, ]) > 0)
        stop("can't combine chains with unequal start/thin")
    attr(x2, "mcpar") <- c(mcpars[1, 1], sum(mcpars[2, ]), mcpars[3,
        1])
    x2
}

## replace names ...
shorten_ag_names <- function(x) {
  rfun <- function(s) {
    ## parental care / sperm competition
    (s %>% stringr::str_replace_all("care", "pc")
      %>% stringr::str_replace_all("spawning", "sc")
    )
  }
  if (is.matrix(x)) {
    colnames(x) <- rfun(colnames(x))
  } else {
    names(x) <- rfun(names(x))
  }
  return(x)
}

## set up contrast matrix - ALL ZEROS, don't overwrite existing file!
setup_contrasts <- function(parnames,
                            cat1 = c("intercept", "pc", "sc", "pcxsc"),
                            cat2 = c("loss", "gain", "netgain")) {
  cnames <- c(outer(cat1, cat2, paste, sep = "_"))
  m <- matrix(0, nrow = length(parnames), ncol = length(cnames),
              dimnames = list(par = parnames, effect = cnames))
  return(m)
}


##' assume compdata is a list with $phylo element
scale_phylo <- function(phy, type = "sumbranches", scale_val = 1) {
  if (type != "sumbranches") stop("only scaling by total sum of branch lengths")
  orig_sumbranches <- sum(phy$edge.length)
  phy$edge.length <- phy$edge.length / orig_sumbranches * scale_val
  attr(phy, "orig_sumbranches") <- orig_sumbranches
  return(phy)
}

#' @param m the 'map' component of a simmap object
get_state_occ_prop <- function(m) {
  s <- unlist(m)
  tapply(s, list(names(s)), sum)
}

image_plot <- function(m, xlab="", ylab = "", sub = "") {
    require("Matrix")
    image(Matrix(m),
          scales=list(x=list(at=seq(ncol(m)),labels=colnames(m),
                             rot=90),
                      y=list(at=seq(nrow(m)),labels=rownames(m))),
          xlab = xlab, ylab = ylab, sub = sub)
}


## how do we incorporate these (elegantly/easily) into the contrasts calculation?

tidy.mle2 <- function (x, conf.int = FALSE, conf.level = 0.95,
                       conf.method = "spline", ...)
{
    co <- bbmle::coef(bbmle::summary(x))
    ret <- broom:::as_tidy_tibble(co, new_names = c("estimate",
        "std.error", "statistic", "p.value"))
    if (conf.int) {
        ci <- bbmle::confint(x, level = conf.level, method = conf.method)
        if (is.null(dim(ci))) {
            ci <- matrix(ci, nrow = 1)
        }
        colnames(ci) <- c("conf.low", "conf.high")
        ci <- dplyr::as_tibble(ci)
        ret <- dplyr::bind_cols(ret, ci)
    }
    dplyr::as_tibble(ret)
}

##' save an MCMC diagnostic pairs plot to a file
##' @param mcmc_obj an mcmc object
##' @param mc_theme ggplot theme (default is black/white with no between-panel spacing
##' @param fn file name for saved plot
##' @param \dots additional arguments to \code{ggsave}
## FIXME: precompute/cache densities?
##' @example
##' library(ggplot2); library(coda)
##' m <- replicate(4, as.mcmc(matrix(rnorm(4000), ncol = 4)), simplify = FALSE)
##' mk_mcmcpairsplot(m, return_val = "plot")
mk_mcmcpairsplot <- function(mcmc_obj, fn, mc_theme = NULL,
                             return_val = c("fn", "plot"), ...) {
  return_val <- match.arg(return_val)
  if (is.null(mc_theme)) {
    mc_theme <- theme_bw() + theme(panel.spacing = grid::unit(0, "lines"))
  }
  theme_set(mc_theme)
  p <- GGally::ggpairs(as.data.frame(emdbook::lump.mcmc.list(mcmc_obj)), progress=FALSE,
               lower=list(continuous=function(...) my_mcmc(..., show_prior=FALSE)),
               upper=list(continuous=function(...) my_mcmc(..., geom="density",
                                                           show_prior = FALSE)))
  if (return_val == "fn") {
      ggsave(filename = fn, plot = p, ...)
      return(fn)
  }
  return(p)
}

fix_cnms <- function(x) {
    dplyr::rename(x, lwr = "conf.low", upr = "conf.high")
}

## generic tidying function/wrapper, allows for incorporating contrasts
my_tidy <- function(x, contrast_mat = NULL, conf.level = 0.95) {
    require("broom")
    require("broom.mixed")
    if (inherits(x, "mcmc.list")) x <- emdbook::lump.mcmc.list(x)
    ## shouldn't have to do this (use method dispatch), but we want
    ##  different arguments for different types:
    if (inherits(x, "mcmc")) {
        if (!is.null(contrast_mat)) x <- coda::as.mcmc(x %*% contrast_mat)
        res <- tidy(x, conf.int = TRUE, robust = TRUE)
    } else if (inherits(x, "corhmm")) {
        res <- tidy(x, conf.int = TRUE, contrast_mat = contrast_mat)
    } else if (inherits(x, "mle2")) {
        if (is.null(contrast_mat)) {
            res <- tidy(x, conf.int = TRUE, conf.method = "quad")
        } else {
            p <- drop(coef(x) %*% contrast_mat)
            V <- t(contrast_mat) %*% vcov(x) %*% contrast_mat
            ## copied from tidy.corhmm() above
            sds <- sqrt(diag(V))
            qq <- qnorm((1+conf.level)/2)
            res <- dplyr::tibble(
                              term = names(p),
                              estimate = p,
                              std.error = sds,
                              statistic = p/sds,
                              conf.low = estimate - qq*std.error,
                              conf.high = estimate + qq*std.error)
        }
    }
    return(fix_cnms(res))
}

## beginnings of a function to do a 'simple' static plot
my_tar_network <- function() {
    require(igraph)
    tt <- tar_network()
    gg <- graph_from_data_frame(tt$edges)
    gg <- topo_sort(gg)
}

gainloss_ind_prs <- function(x) {
    imat <- x$index.mat
    upr <- unique(na.omit(imat[upper.tri(imat)]))
    lwr <- unique(na.omit(imat[lower.tri(imat)]))
    Map(c, upr, lwr)
}

profile.corhmm <- function(fitted, max_val = 3, delta = 0.1, maxit = 50,
                           params = NULL,
                           verbose = FALSE, optControl = list(maxit = 20000),
                           conv_action = stop, ...) {
    ## FIXME: allow parameter selection?
    ## *generic* profiling recipe?
    f <- make_nllfun(fitted)
    p0 <- p <- fitted$args.list$p
    L0 <- -1*c(logLik(fitted))
    ## FIXME: check for sd problems?
    sd <- sqrt(diag(solve(attr(fitted, "hessian"))))
    offset <- sd*delta
    res <- list()
    ## can't use 'par' as loop variable, conflicts with optim() argument name
    if (is.null(params)) params <- seq_along(p0)
    for (cur_par in params) {
        if (verbose) cat(sprintf("param %d (%s)\n", cur_par, names(p)[cur_par]))
        ## construct full parameter vector from restricted parameters
        mkpar <- function(p, cur_par, cur_parval) {
            pp <- rep(NA, length(p)+1)
            pp[-cur_par] <- p
            pp[cur_par] <- cur_parval
            names(pp) <- names(p0)
            pp
        }
        ## evaluate negative log-lik over restricted parameters
        wrapfun <- function(p, cur_par, cur_parval) {
            f(mkpar(p, cur_par, cur_parval))
        }
        ## FIXME: parallel evaluation?
        for (dir in c(-1, 1)) {
            if (verbose) cat("dir: ", dir, "\n")
            cur_parval <- p0[cur_par]
            p <- p0
            it <- 0
            val <- 0
            while (it < maxit && val < max_val) {
                cur_parval  <- cur_parval + dir*offset[cur_par]
                if (verbose) cat(cur_parval, " ")
                opt <- optim(fn = wrapfun, par = p[-cur_par], cur_par = cur_par, cur_parval = cur_parval,
                             control = optControl)
                if (opt$convergence != 0) conv_action("convergence code ", opt$convergence)
                val <- opt$value - L0
                p <- mkpar(p, cur_par, cur_parval)
                if (verbose) cat(val, "\n")
                res <- c(res, list(
                                  data.frame(
                                      .zeta = -2*sqrt(val)*dir,
                                      .par = names(p0)[[cur_par]],
                                      .focal = cur_parval,
                                      rbind(p)
                                  )  ## data.frame()
                              )  ## list()
                         ) ## c()
                it <- it + 1
            }  ## while goal not achieved
        }  ## loop over dir
    } ## loop over params
    ##  assemble results
    do.call("rbind", res)
}
## debug(profile.corhmm)
## profile(fitted, verbose = TRUE)

fit_contrast.corhmm <- function(fitted, contrast, fixed_vals,
                                optControl = list(maxit = 20000),
                                ## HACK: contrasts won't be on quite the same scale as
                                ##  the original transitions ...
                                lower = log(0.001),
                                upper = log(10000)) {
    require("bbmle")
    f <- make_nllfun(fitted)
    p0 <- coef(fitted)
    ## updated contrast matrix (transposed & reordered)
    c2 <- t(contrast)[,names(p0)]
    ## reference (starting) values in contrast space
    p0_c <- drop(c2 %*% p0)
    ## inverse-contrast matrix 
    ic <- solve(c2)
    stopifnot(all.equal(p0, drop(ic %*%  p0_c)))
    ## DRY?
    mkpar <- function(p, cur_par, cur_parval) {
        pp <- rep(NA, length(p)+length(cur_par))
        pp[-cur_par] <- p
        pp[cur_par] <- cur_parval
        names(pp) <- names(p0)
        pp
    }
    ## evaluate negative log-lik over restricted parameters
    wrapfun <- function(p, cur_par, cur_parval) {
        pp <- mkpar(p, cur_par, cur_parval)
        f(ic %*% pp)
    }
    cur_par <- match(names(fixed_vals), names(p0_c))
    if (length(cur_par) == 0) stop("bad param names")
    start <- p0_c[-cur_par]
    parnames(wrapfun) <- names(start)
    ## wrapfun(start, cur_par = cur_par, cur_parval = fixed_vals)
    ## maybe not necessary/env?
    w2 <- function(p) wrapfun(p, cur_par = cur_par, cur_parval = fixed_vals)
    fit <- nloptwrap(par = start,
                     fn = w2,
                     lower = lower,
                     upper = upper)
    nm <- setdiff(colnames(contrast), names(fixed_vals))
    names(fit$par) <- nm
    ## HACK: create an mle2 object (mle2 doesn't take user-specified optimizers,
    ##  and I want to use nloptwrap ...)
    H <- optimHess(fit$par, fn = w2)
    dimnames(H) <- list(nm, nm)
    ret <- c(fit, list(hessian = H))
    class(ret) <- c("corhmm_contrast", "list")
    return(ret)
}

logLik.corhmm_contrast <- function(x) {
    LL <- -1*x$fval
    df <- length(x$par)
    attr(LL, "df") <- df
    class(LL) <- "logLik"
    return(LL)
}

tidy.corhmm_contrast <- function(x,
                        conf.int = FALSE,
                        conf.method = "wald",
                        conf.level = 0.95,
                        p.value = FALSE,
                        profile = NULL,
                        exponentiate = FALSE,
                        prof_args = NULL,
                        contrast_mat = NULL,
                        ...) {
    p <- x$par
    res <- dplyr::tibble(term = names(p), estimate = p)
    if (conf.int) {
        V <- solve(x$hessian)
        sds <- sqrt(diag(V))
        qq <- qnorm((1+conf.level)/2)
        res <- dplyr::mutate(res,
                             std.error = sds,
                             statistic = p/sds,
                             conf.low = estimate - qq*std.error,
                             conf.high = estimate + qq*std.error)
    }
    if (exponentiate) {
        res <- mutate(res, across(estimate, exp))
        res <- mutate(res, across(starts_with("conf"), exp))
        if ("std.error" %in% names(res)) {
            res <- mutate(res, across(std.error, ~ . * estimate))
        }
    }
    return(res)
}

nloptwrap <- function (par, fn, lower = -Inf, upper = Inf, control = list(), ...) {
    defaultControl <- list(algorithm = "NLOPT_LN_BOBYQA", xtol_abs = 1e-08, ftol_abs = 1e-08, 
                           maxeval = 1e+05)
    lower <- rep(lower, length.out = length(par))
    upper <- rep(upper, length.out = length(par))
    for (n in names(defaultControl)) if (is.null(control[[n]])) 
                                         control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0 = par, eval_f = fn, lb = lower, ub = upper, 
                  opts = control, ...)
    with(res, list(par = solution, fval = objective, feval = iterations, 
                   conv = if (status < 0 || status == 5) status else 0, 
                   message = message))
}

if (FALSE) {
    library(corHMM)
    library(targets)
    tar_load(ag_model_pcsc)
    source("R/utils.R")
    ## read contrast matrix and convert to matrix
    invcontr <- read.csv("contr_invertible.csv")
    rn <- invcontr$parname
    invcontr <- as.matrix(invcontr[,-1])
    dimnames(invcontr) <- list(rn, colnames(invcontr))
    ## fit additive model
    fit_contrast.corhmm(ag_model_pcsc, invcontr,
                        ## zero out interactions
                        fixed_vals = c(pcxsc_loss = 0, pcxsc_gain = 0),
                        optControl = list(maxit = 20000, trace = 2))
}

latex_mat <- function(m, align = "r") {
    hdr <- c("\\left[",
             sprintf("\\begin{array}{%s}", paste(rep(align, ncol(m)), collapse = "")))
    body <- character(nrow(m))
    for (i in seq(nrow(m))) {
        body[i] <- paste(paste(m[i,], collapse = " & "), "\\\\")
    }
    tail <- c("\\end{array}", "\\right]")
    c(hdr, body, tail)
}

mk_flowfig <- function(model = ag_model_pcsc, tikz = FALSE,
                       nmag = 0.04, with_labs = FALSE) {
    M <- model$solution
    R <- model$args.list$rate
    dimnames(R) <- dimnames(M)
    col_vec1 <- rep(NA, max(R)-1)
    R[R==max(R)] <- NA

    ## constrained model colour vectors
    ## 12-parameter model
    ## match by name ...
    ## up_down_ind == all arrows between equal ag*
    ag_labs <- substr(colnames(R), 1, 3)
    updown <- outer(ag_labs, ag_labs, "==")
    up_down_ind <- sort(unique(R[!is.na(R) & updown])) ## c(1,2,4,6) 
    col_vec1[up_down_ind] <- gray(c(0, 0.3, 0.6, 0.9))
    col_vec1[-up_down_ind] <- qualitative_hcl(8)
    C1 <- c(R)
    C1 <- col_vec1[C1] ## map colors to indices
    dim(C1) <- dim(R)

    ## pc-only model
    col_vec2 <- col_vec1
    pc_pairs <- list(c(5, 6), c(7, 8), c(1, 2), c(3, 4))
    for (i in seq_along(pc_pairs)) {
        col_vec2[pc_pairs[[i]][2]] <- col_vec2[pc_pairs[[i]][1]]
    }
    C2 <- c(R)
    C2 <- col_vec2[C2] ## map colors to indices
    dim(C2) <- dim(R)

    ## sc-only model
    col_vec3 <- col_vec1
    sc_pairs <- list(c(5, 7), c(6, 8), c(1, 3), c(2, 4))
    for (i in seq_along(sc_pairs)) {
        col_vec3[sc_pairs[[i]][2]] <- col_vec3[sc_pairs[[i]][1]]
    }
    C3 <- c(R)
    C3 <- col_vec2[C3] ## map colors to indices
    dim(C3) <- dim(R)

    mksymm <- function(X) {
        eqzero <- function(z) !is.na(z) && z == 0
        for (i in 1:nrow(X)) {
            for (j in 1:ncol(X)) {
                if (eqzero(X[i,j]) && !eqzero(X[j,i])) {
                    X[i,j] <- X[j,i]
                }
                if (eqzero(X[j,i]) && !eqzero(X[i,j])) {
                    X[j,i] <- X[i,j]
                }
            }
        }
        X
    }

    ## nudges for arrows
    Nx <- M
    Nx[!is.na(Nx)] <- 0
    Ny <- Nx

    ## outer up/down arrows: nudge L/R
    Nx[matrix(c("ag1_pc1_sc0", "ag1_pc0_sc0",
                "ag1_pc1_sc1", "ag1_pc0_sc1"),
              byrow = TRUE, ncol = 2)] <- c(-nmag, nmag)
    Nx <- mksymm(Nx)
    ## outer L/R arrows: nudge up/down
    Ny[matrix(c("ag1_pc1_sc0", "ag1_pc1_sc1",
                "ag1_pc0_sc0", "ag1_pc0_sc1"),
              byrow = TRUE, ncol = 2)] <- c(nmag, -nmag)
    Ny <- mksymm(Ny)

    ## turn off arrow labels
    storage.mode(R) <- "character"
    R[!is.na(R)] <- ""

    ## node arrangement
    ## hack (reverse order)
    vals <- expand.grid(sc=0:1, pc=0:1, ag=0:1)[,3:1]
    ## put nodes on flattened hypercube
    ##  (ag=0 inner square, ag = 1 outer square

    box_mag <- vals$ag+0.5
    shift <- c(0.3, 0.5)
    xval <- with(vals, (sc-0.5)*box_mag)
    yval <- with(vals, (pc-0.5)*box_mag)
    pos <- cbind(xval, yval)
    pos2 <- (pos-min(pos) + shift[1])/(diff(range(pos))*(1+shift[2]))

    
    ## labels on subset of arrows (for constraint example)
    R2 <- R
    from <- sprintf("ag0_pc%d_sc%d", rep(0:1,2), rep(0:1, each=2))
    to <-   sprintf("ag1_pc%d_sc%d", rep(0:1,2), rep(0:1, each=2))
    if (!tikz) {
        R2[cbind(from, to)] <- LETTERS[1:4]
    } else {
        R2[cbind(from, to)] <- sprintf("\\encircle{%s}", LETTERS[1:4])
    }

    ## protect
    ## (don't need if latex=TRUE)
    if (tikz) {
         rownames(R) <- rownames(R2) <- gsub("_", ".", rownames(R))
         colnames(R) <- colnames(R2) <- gsub("_", ".", colnames(R))
    }

    mkplot <- function(mat = R, C = C1) {
        plotmat(mat, pos = pos2, xlim = c(-3,3),
                ## arr.lwd = sqrt(M/50),
                box.type = "ellipse", box.prop = 0.5,
                arr.lcol = C, arr.col = C,
                arr.nudge.x = Nx,
                arr.nudge.y = Ny,
                latex = tikz)
    }
    
    if (!with_labs) {
        mkplot()
    } else {
        ## with labels, for supp/contrast explanation
        pp <- mkplot(R2, C1)
        w <- nzchar(na.omit(c(R2)))
        if (!tikz) {
            ## not quite sure why this fussing is required
            nx <- 0.015*c(0.75,0.75,-0.75,-1)
            ny <- 0.015*c(-0.9,0.9,0.9,0)
            with(pp$arr, mapply(plotrix::draw.circle, TextX[w]+nx, TextY[w]+ny,
                                MoreArgs=list(radius=0.025)))
        }
    }
    invisible(NULL)
}

tikz <- function(file, ...) {
    tikzDevice::tikz(file, ..., standAlone = TRUE)
}

enc_line <- "\\newcommand\\encircle[1]{ \\tikz[baseline=(X.base)]  \\node (X) [draw, shape=circle, inner sep=0] {\\strut #1};}"

## need tikzDevice to create these definitions, but want to define pkg_install() ...
if (require("tikzDevice")) {
    tLP <- getOption("tikzLatexPackages")
    if (!any(grepl("encircle", tLP))) {
        options("tikzLatexPackages" = c(tLP,
                                        enc_line))
    }
}

## tikz_insert <- function(file, str = enc_line, line = 12) {
##     system(sprintf("sed '%d i %s' %s", line = , str, file))
## }


get_chains <- function(results, ret_val = c("data.frame", "mcmc", "mcmc.list"), xcols = c("chain", "Iteration"),
                       maxlen = 6000) {
    ret_val <- match.arg(ret_val)
    rate_cols <- grep("^q[0-9]+", colnames(results$Log$results), value = TRUE)
    chains <- results$Log$results[,c(xcols, rate_cols)]  ## q** values only
    if ((nc <- nrow(chains)) > maxlen) {
        message(sprintf("keeping %d of %d rows", maxlen, nc))
        keep <- round(seq(1, nc, length.out = maxlen))
        chains <- chains[keep,]
    }
    cols_disallowed <- which(apply(chains==0, 2, all)) ## forbidden/boring
    dupes <- c("q24","q57","q68", ## care gain
               "q42","q75","q86", ## care loss
               "q34","q56","q78", ## spawn gain
               "q43","q65","q87" ## spawn loss
               )
    cols_dupes <- match(dupes, colnames(chains))
    res <- chains[, -c(cols_dupes, cols_disallowed)]
    res <- switch(ret_val,
                  mcmc = as.mcmc(res[!names(res) %in% xcols]),
                  mcmc.list = as.mcmc.list(lapply(split(res[!names(res) %in% xcols], res$chain), mcmc)),
                  res)
    return(res)
}

## computing contrasts - from rates to geometric mean ratios
gmean <- function(x, y) sqrt(x*y)  ## equivalent to exp((log(x) + log(y))/2)
cfun_nonlog <- function(q1, q2, q3, q4) {
    gmean(q1,q2)/gmean(q3,q4)
    ## or exp(
    ##        (log(q1) + log(q2))/2 -
    ##        (log(q3) + log(q4))/2
    ##    ) 
}

## function for getting contrasts from rates
get_contrasts <- function(results, maxlen = 6000) {
    chains <- results$Log$results
    ## DRY, this is also in get_chains()
    if ((nc <- nrow(chains)) > maxlen) {
        message(sprintf("keeping %d of %d rows", maxlen, nc))
        keep <- round(seq(1, nc, length.out = maxlen))
        chains <- chains[keep,]
    }
    chains |> 
        mutate(gain_care_effect =   cfun_nonlog(q37,q48,q15,q26),
               gain_spawn_effect =  cfun_nonlog(q26,q48,q15,q37),
               gain_interaction =   cfun_nonlog(q15,q48,q37,q26),
               loss_care_effect =   cfun_nonlog(q73,q84,q51,q62),
               loss_spawn_effect =  cfun_nonlog(q62,q84,q51,q73),
               loss_interaction =   cfun_nonlog(q51,q84,q73,q62))
}

## should be covered by BMB fork of bayestestR now
## ## shims for bayestestR
## abind <- function(x) {
##   nchains <- length(x)
##   a <- array(dim = c(nrow(x[[1]]), nchains, ncol(x[[1]])))
##   for (i in seq(nchains)) {
##     a[,i,] <- as.matrix(x[[i]])
##   }
##   if (!is.null(parnames <- dimnames(x[[1]])[[2]])) {
##     dimnames(a)[[3]] <- parnames
##   }
##   return(a)
## }

## #' @examples
## #' mkdata <- function(nrow = 1000, ncol = 2, parnm = LETTERS[1:ncol]) {
## #'    replicate(ncol, rnorm(nrow) |> as.data.frame() |> setNames(parnm)
## #' }
## #' dd <- replicate(5, mkdata(), simplify = FALSE)
## #' x <- diagnostic_posterior(dd)
## diagnostic_posterior.default <- function(posterior, ...) {
##   if (is.list(posterior)) {
##     for (i in seq.along(posterior)) {
##       p <- posterior[[i]]
##       if (!(
##         (inherits(p, "data.frame") || inherits(p, "matrix")) &&
##         (dim(p) == 2) &&
##         (ncol(p) == ncol(posterior[[1]])))
##         ) stop("'posterior' must be a 3D array or a list of data frames with equal numbers of columns")
##     }
##     posterior <- abind(posterior)
##   }
##   if (!(inherits(posterior, "array") && length(dim(posterior)) == 3)) {
##     stop("expecting a 3D array")
##   }
##   ret <- rstan::monitor(posterior)
##   class(ret) <- "data.frame" ## get rid of simsummary() class
##   ret <- ret |> tibble::rownames_to_column(var = "Parameter")
##   as_tibble(ret)
## }
