#### BayesTraits Model ####

## should modularize this some more:
##   allow path/bin/cores to be set externally
##   split into individual `make` or `targets` targets so they
##    can be (re)run separately if necessary

## FIXME/TODO:
##  * finish testing/considering RJ results
##     * compare posterior dist of rates (rather than contrasts),
##       to avoid dealing with log-scale contrasts involving zero?
##  * simmap? time-weighted contrasts?
##  * Compare BT 'regular' with our results/side-by-side plots
##  * more Bayes diagnostics - improved R-hat? rank-histogram plots from
##     bayestestR?

## command vector
## use '-Threaded' for threaded executable, specify number of cores
options(bt_path = "BayesTraitsV4.0.0-Linux-Threaded",
        bt_bin = "BayesTraitsV4")

options(bayestraits.cores = 16)
## run from head directory of repo
## (for target loading, readRDS)

## packages
library(btw)
library(targets)
library(ape)
library(dplyr)

## needed for testing only

library(ggplot2); theme_set(theme_bw())
library(tidyr)

source("R/utils.R") ## for get_chains

## loading trees
tar_load(treeblock)
trees <- do.call(c, treeblock)
trees <- .compressTipLabel(trees)

## loading data

tar_load(ag_compdata_tb)

## this is clever and extensible but longer than the brute-force original ...
codefun <- function(ag, pc, sc) {
    ## octal coding (more or less) by position
    v <- c(ag, pc, sc)
    f <- function(z) sum(c(4,2,1)*as.numeric(z)) + 1
    nq <- sum(v=="?")
    ## replace *first* ? with i
    ssub <- function(i, x) {
        x[which(x=="?")[1]] <- i
        return(x)
    }
    if (nq ==0) return(f(v))
    if (nq ==1) {
        return(as.numeric(
            paste0(f(ssub(0, v)),
                   f(ssub(1,v)))
        ))
    }
    ## handle nq==2 case cleverly (mapping over 0/1 combos) or by brute force?
    ## could use sub() to substitute 0/1 for *first* ?, then *second* ? ...
    fsub <- function(i,j, x) {
        f(ssub(j, ssub(i, x)))
    }
    return(as.numeric(
        paste0(fsub(0,0,v),
               fsub(0,1,v),
               fsub(1,0,v),
               fsub(1,1,v))
    ))
}

## testing
## codefun(0, "?", 0)
## codefun(0, "?", "?")
    
tdata <- ag_compdata_tb$data %>%
    mutate(state=factor(case_when(
    (ag==0 & pc==0 & sc==0) ~ 1,
    (ag==0 & pc==0 & sc==1) ~ 2,
    (ag==0 & pc==1 & sc==0) ~ 3,
    (ag==0 & pc==1 & sc==1) ~ 4,
    (ag==1 & pc==0 & sc==0) ~ 5,
    (ag==1 & pc==0 & sc==1) ~ 6,
    (ag==1 & pc==1 & sc==0) ~ 7,
    (ag==1 & pc==1 & sc==1) ~ 8,
    (ag==0 & pc==0 & sc=="?") ~ 12,
    (ag==0 & pc==1 & sc=="?") ~ 34,
    (ag==1 & pc==0 & sc=="?") ~ 56,
    (ag==1 & pc==1 & sc=="?") ~ 78,
    (ag==0 & pc=="?" & sc==0) ~ 13,
    (ag==0 & pc=="?" & sc==1) ~ 24,
    (ag==1 & pc=="?" & sc==0) ~ 57,
    (ag==1 & pc=="?" & sc==1) ~ 68,
    (ag==0 & pc=="?" & sc=="?") ~ 1234,
    (ag==1 & pc=="?" & sc=="?") ~ 5678))) %>% 
    dplyr::select(species,state)

summary(tdata)

#### PRIOR SCALING ####

## where is this from??? presumably described/justified in original
##  writeup?
rate_prior_0 <- rate_prior <- c(4.236, 1.41)
## BT scales branch lengths to a *mean* of 0.1 by default

edge_len_tot <- sapply(trees, \(x) sum(x$edge.length))
stopifnot(all(abs(edge_len_tot - 1) < 1e-8))
n_edges <- sapply(trees, \(x) nrow(x$edge))  ## all 1212
stopifnot(var(n_edges) == 0)
n_edge <- n_edges[[1]]

prior_fac <- (1/n_edge)/0.1
## so mean branch length in our analysis is 1/1212
## our branch lengths are much *shorter*
## therefore our rates will be *higher* than BT's
## therefore we should decrease the priors we give to BT by log(prior_fac)
## (-4.79)
rate_prior[1] <- exp(rate_prior_0[1] + log(prior_fac))
prior1 <- sprintf("PriorAll lognormal %f %f", rate_prior[1], rate_prior[2])
prior2 <- sprintf("RevJump lognormal %f %f", rate_prior[1], rate_prior[2])
prior3 <- sprintf("RevJump exp 10") ## match 'traditional' RJ prior

#### COMMAND FUNCTION ####

##  q14 q16 q17 q18 q23 q25 q27 q28 q32 q35 q36 q38 q41 q45 q46 q47 q52 q53 q54 q58 q61 q63 q64 q67 q71 q72 q74 q76 q81 q82 q83 q85 0", ## impossible rates 
zero_rates <- c(14, 16:18, 23, 25, 27:28, 32, 35:36, 38, 41, 45:47, 52:54,
                58, 61, 63:64, 67, 71:72, 74, 76, 81:83, 85)
qz <- paste("q", zero_rates, sep ="", collapse = " ")

bt_command <- function(prior = NULL,
                       iterations = 51e4,
                       burnin = 1e4,
                       sample = 1,
                       seed = 101,
                       cores = getOption("bayestraits.cores", 1)) {
    cvec <- c("1", ## MultiState
             "2", ## MCMC
              "ScaleTrees", ## scaling branch lengths to a mean of 0.1
              "AddTag Root Erpetoichthys_calabaricus Mugil_liza", ## adding a tag at the root
              "Fossil Node01 Root 2", ## fossilizing the root (fix trait value at root: group spawning, no pc/ag)
              paste("Res", qz, "0"),
              "Res q13 q24 q57 q68", ## care gain
              "Res q31 q42 q75 q86", ## care loss
              "Res q12 q34 q56 q78", ## spawn gain
              "Res q21 q43 q65 q87"  ## spawn loss
              )
    if (!is.null(prior)) cvec <- c(cvec, prior)
    cvec <-  c(cvec,
               sprintf("Iterations %d", iterations),
               sprintf("Burnin %d", burnin),
               sprintf("Seed %d", seed),
               sprintf("Sample %d", sample)  ## thinning
               )
    if (!is.null(cores)) cvec <- c(cvec, sprintf("Cores %d", cores))
    return(cvec)
}

## running in parallel via parLapply etc. would be tricky because of hard-coded output files -- would
##  need to figure out how to get btw::bayestraits() to disambiguate/send them different places
bt_run <- function(data = NULL, trees = NULL, prior = NULL, dir = "bayestraits",
                   fn = "", chains = 4, seed0 = 100, verbose = FALSE,
                   ...) {
    all_res <- lapply(seq(chains),
                      function(i) {
                          if (verbose) cat("chain ", i, "\n")
                          command <- bt_command(prior = prior, ...,
                                                seed = seed0 + i)
                          tt <- system.time(
                              res <- bayestraits(data, trees, command)
                          )
                          attr(res, "time") <- tt
                          return(res)
                      })
    results <- list(Log = list(options = all_res[[1]]$Log$options, ## same opts for all chains
                               results = do.call(rbind, lapply(seq_along(all_res), function(i) data.frame(chain = i, all_res[[i]]$Log$results)))
                               ))
    results[["time"]] <- lapply(all_res, function(x) attr(x, "time"))
    for (c in c("Schedule", "Stones", "AncStates", "OutputTrees")) {
        results[[c]] <- lapply(all_res, function(x) x[[c]])
    }
    if (nchar(fn) > 0) saveRDS(results, file = file.path(dir, fn))
    return(results)
}

#### DATA, REGULAR, DEFAULT PRIORS ####

r1 <- bt_run(data = tdata, trees, fn = "bt_model_data_reg_default.rds", verbose = TRUE)

## check
r1B <- readRDS("bayestraits/bt_model_data_reg_default.rds")
r1B_L <- (r1B
    |> get_chains()
    |> tidyr::pivot_longer(-c(chain, Iteration))
)
ggplot(r1B_L, aes(Iteration, value, colour = factor(chain))) + geom_line() + facet_wrap(~name, scale = "free") +
    scale_y_log10()

##### TDATA, REGULAR, OUR PRIORS ####

r2 <- bt_run(prior = prior1, data = tdata, trees, fn = "bt_model_data_reg_priors.rds", verbose = TRUE)

#### DATA-LESS, REGULAR, DEFAULT PRIORS ####

## data-less data
data_dataless <- mutate(tdata, state=12345678)
summary(data_dataless)

r3 <- bt_run(prior = NULL, data = data_dataless, trees, fn = "bt_model_nodata_reg_default.rds", verbose = TRUE)

#### DATA-LESS, REGULAR, OUR PRIORS ####

r4 <- bt_run(prior = prior1, data = data_dataless, trees, fn = "bt_model_nodata_reg_priors.rds", verbose = TRUE)

#### TDATA, RJ, DEFAULT ####

r5 <- bt_run(prior = "revJump uniform 0 100",
             data = tdata, trees, fn = "bt_model_data_rj_default.rds", verbose = TRUE)

#### DATA-LESS, RJ, DEFAULT ####

r6 <- bt_run(prior = "revJump uniform 0 100",
             data = data_dataless, trees, fn = "bt_model_nodata_rj_default.rds", verbose = TRUE)

#### TDATA, RJ, OUR PRIORS ####

r7 <- bt_run(prior = prior2,
             data = tdata, trees, fn = "bt_model_data_rj_priors.rds", verbose = TRUE)

#### DATA-LESS, RJ, OUR PRIORS ####

r8 <- bt_run(prior = prior2,
             data = data_dataless, trees, fn = "bt_model_nodata_rj_priors.rds", verbose = TRUE)

## ALTERNATIVE RJ prior (exp(10))

## based on trace plots, seems just as bad as the other RJ+data options
## diagnostic plots fail because ESS calculation fails (Chol inversion)
if (FALSE) {
    r9 <- bt_run(prior = prior3,
                 data = tdata, trees, fn = "bt_model_data_rj_prior-exp10.rds", verbose = TRUE)

    ## looking better?
    r10 <- bt_run(prior = prior3,
                  chains = 1,
                  data = tdata,
                  trees,
                  fn = "bt_model_data_rj_prior-exp10-long.rds",
                  verbose = TRUE,
                  ## passed to bt_command:
                  iterations = 1e7,
                  sample = 100
                  )

    r11 <- bt_run(prior = prior2,
                  chains = 4,
                  data = tdata,
                  trees,
                  fn = "bt_model_data_rj_prior2-long.rds",
                  verbose = TRUE,
                  ## passed to bt_command:
                  iterations = 1e7,
                  sample = 100
                  )
    
    r12 <- bt_run(prior = prior2,
                  chains = 4,
                  data = data_dataless,
                  trees,
                  fn = "bt_model_nodata_rj_prior2-long.rds",
                  verbose = TRUE,
                  ## passed to bt_command:
                  iterations = 1e7,
                  sample = 100
                  )

}
