# loading packages
if (packageVersion("corHMM") < "2.6.1") {
  stop("need hacked version of corHMM. Try 'remotes::install_github(\"bbolker/corHMM\", build_vignettes=TRUE)'")
}
## if installation fails when building vignettes drop build_vignettes=TRUE (default is FALSE)
library(corHMM)
library(bbmle) ## for likelihood profiles/CIs
library(numDeriv) ## for Hessian/Wald CIs
library(MCMCpack) ## for Bayes (Metropolis-Hastings) run
library(coda)
library(lattice)
library(ggplot2)
theme_set(theme_bw())
library(colorspace)
library(GGally)
library(ggthemes)
library(fishtree)
library(caper)
### remotes::install_github("bbolker/corHMM")

# loading data and tree and trimming
allData <- read.csv("binaryTraitData.csv", header = TRUE)
fullPhy <- fishtree_phylogeny(allData$species)
traitData <- data.frame(species = allData$species, ag = allData$ag, care = allData$care, mating = allData$mating, names = allData$species)
trimmedData <- comparative.data(fullPhy, traitData, names.col = names)
phy <- trimmedData$phy
data <- trimmedData$data

## default model
if (FALSE) {
  (MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1))
}

# getting state matrix and legend
LegendAndRateMat <- getStateMat4Dat(data)
RateMat <- LegendAndRateMat$rate.mat
getStateMat4Dat(data)
LegendAndRateMat <- getStateMat4Dat(data)
RateMat <- LegendAndRateMat$rate.mat

# making the gain and loss of care and sperm comp exogenous, down to 12 pars
pars2equal <- list(c(7, 10, 20, 23), c(4, 11, 17, 24), c(2, 5, 15, 18), c(1, 8, 14, 21))
StateMatA_constrained <- equateStateMatPars(RateMat, pars2equal)
StateMatA_constrained

# making the simplified model
MK_3state_simple <- corHMM(phy = phy, data = data, rate.cat = 1, rate.mat = StateMatA_constrained)
MK_3state_simple

# Check that we can retrieve the same log-likelihood as above.
do.call(corHMM:::dev.corhmm, MK_3state_simple$args.list)

# Make a wrapper function that calls dev.corhmm with a specified parameter vector, and all the other stuff it needs:
nllfun <- function(p) {
  a <- MK_3state_simple$args.list
  a$p <- p
  do.call(corHMM:::dev.corhmm, a)
}

# Compute Wald confidence intervals (standard-error based)
p <- MK_3state_simple$args.list$p
H <- numDeriv::hessian(nllfun, p)
sds <- sqrt(diag(solve(H)))
wald.ci <- sweep(qnorm(0.975) * outer(sds, c(-1, 1)), 1, FUN = "+", STATS = p)

names(p) <- parnames(nllfun) <- paste("p",1:12)
nllfun(p)

# fit an mle2 model
m0 <- mle2(minuslogl=nllfun, start=p, vecpar=TRUE, lower=log(1e-9), upper=log(100), method="L-BFGS-B")

print(system.time(pp <- profile(m0,std.err=0.2,trace=TRUE)))

plot(pp, show.points = TRUE)

prof.ci <- confint(pp)

# Bayesian
logpostfun <- function(p, lb = log(1e-9), ub = log(1e2), range = 3) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  loglik <- -1 * nllfun(p)
  log.prior <- sum(dnorm(p, mean = prior.mean, sd = prior.sd, log = TRUE))
  return(loglik + log.prior) ## product of likelihood and prior -> sum of LL and log-prior
}

fn <- "MK_3state_simple_mcmc.rds"
if (!file.exists(fn)) {
  m1 <- MCMCmetrop1R(logpostfun, p, verbose = 1000, mcmc = 100000)
  saveRDS(m1, file = fn)
} else {
  m1 <- readRDS(fn)
}

xyplot(m1)
summary(m1)
raftery.diag(m1)
pairs(as.matrix(m1), gap = 0, pch = ".")
bayes.ci <- t(apply(m1, 2, quantile, c(0.025, 0.975)))

# graphing all CIs
all.ci <- setNames(
  as.data.frame(rbind(wald.ci, prof.ci, bayes.ci)),
  c("lwr", "upr")
)

all.ci <- data.frame(
  rate = rep(factor(1:4), 3),
  method = rep(c("Wald", "profile", "Bayes"), each = 4),
  est = c(rep(p, 2), colMeans(m1)),
  all.ci
)
all.ci$upr[is.na(all.ci$upr)] <- Inf ## extend NA confidence interval to limit of graph

ggplot(
  all.ci,
  aes(rate, est, ymin = lwr, ymax = upr, colour = method)
) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_colour_discrete_qualitative() +
  coord_flip()
