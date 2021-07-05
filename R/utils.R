pkgList <- c("tidyverse", "bbmle", "coda", "numDeriv",
             "ggthemes", "fishtree", "caper", "broom.mixed",
             "emdbook", "ramcmc", "corHMM",
             "GGally", "colorspace", "ggmosaic", "targets", "tarchetypes",
             "abind")

## install uninstalled pkgs from pkgList
## check corHMM version, install from bb repo if necessary
## FIXME: save required version elsewhere/externally?
install_pkgs <- function() {
    ip <- installed.packages()
    to_install <- setdiff(pkgList, c(rownames(ip), "corHMM"))
    install.packages(to_install)
    if (!require("corHMM") || packageVersion("corHMM") < "2.7.1.1") {
        remotes::install_github("bbolker/corHMM")
    }
    return(NULL)
}

load_pkgs <- function() {
  invisible(lapply(pkgList, library, character.only = TRUE))
}

#' utility function for hexbin panels for MCMC pairs plot
#' @param data ...
#' @param mapping ...
#' @param lb lower bound
#' @param ub upper bound
#' @param range distance from center to bound, in SD
#' @param show_prior include representation of prior (red lines)
my_mcmc <- function(data, mapping, lb=log(1e-9), ub=log(1e2), range=2,
                    show_prior=FALSE, geom=c("hexbin","density"),
                    density_alpha=0.5,
                    hpd_levels = c(0.5, 0.8, 0.9, 0.95),
                    ...) {
  require("emdbook") ## FIXME: not sure why? remove & see what breaks?
  geom <- match.arg(geom)
  gg1 <- ggplot(data = data, mapping=mapping)
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



image.corhmm <- function(x,
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


## wrap corhmm results to create a negative log-likelihood function
## (faster than calling corHMM(p = exp(log_p))
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
#' @param range_gainloss
#' @param nllfun \emph{negative} log-likelihood function
corhmm_logpostfun <- function(p,
                              lb = log(1e-9),
                              ub = log(1e2),
                              range = 3,
                              gainloss_pairs = NULL,
                              lb_gainloss = log(1e-3),
                              ub_gainloss = log(1e3),
                              range_gainloss = 3,
                              nllfun,
                              treeblock = NULL
                              ) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  loglik <- -1 * nllfun(p)
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
augment_model <- function(model) {
    sn <- state_names(model$data[, -1])
    dimnames(model$solution) <- dimnames(model$index.mat) <- list(sn, sn)
    names(model$args.list$p) <- par_names(model)
    return(model)
}

tidy.corhmm <- function(x,
                        conf.int = FALSE,
                        conf.method = "wald",
                        conf.level = 0.95,
                        p.value = FALSE,
                        profile = NULL,
                        exponentiate = FALSE,
                        prof_args = NULL,
                        ...) {
  f <- make_nllfun(x)
  p <- x$args.list$p

  ## FIXME: augment model here ? or assume already augmented?
  res <- tibble(term = names(p),
                estimate = p)
  if (conf.int) {
    if (conf.method == "wald") {
      H <- numDeriv::hessian(f, p)
      sds <- sqrt(diag(solve(H)))
      qq <- qnorm((1+conf.level)/2)
      res <- mutate(res,
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
  image(Matrix(m),
      scales=list(x=list(at=seq(ncol(m)),labels=colnames(m),
                         rot=90),
                  y=list(at=seq(nrow(m)),labels=rownames(m))),
      xlab = xlab, ylab = ylab, sub = sub)
}


## how do we incorporate these (elegantly/easily) into the contrasts calculation?
