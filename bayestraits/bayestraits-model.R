#### BayesTraits Model analysis ####

## see http://bbolker.github.io/bbmisc/bayes/examples.html
## for more info on convergence diagnostics etc.

## TO DO: improved (Vehtari et al.) R-hat?
## improved trace plots (rank-based) ?

## compare RJ posteriors (for rates, not log-contrasts) with
##  'regular' and ours

## run from head directory of repo
## (for target loading, readRDS)

## packages
library(btw)
library(targets)
library(tidyverse)
library(ape)
library(coda)
library(bayestestR)
library(ggplot2); theme_set(theme_bw())

min_val <- 1e-4 ## for plotting traces etc.: zeros lost

source("R/utils.R")

name_match_vec <- c(loss.sc="q21",
                loss.pc="q31",
                loss.ag_pc0_sc0 = "q51",
                gain.sc = "q12",  ## also q87
                loss.ag_pc0_sc1  = "q62",
                gain.pc = "q13",
                loss.ag_pc1_sc0 = "q73",
                loss.ag_pc1_sc1 = "q84",
                gain.ag_pc0_sc0 = "q15",
                gain.ag_pc0_sc1 = "q26",
                gain.ag_pc1_sc0 = "q37",
                gain.ag_pc1_sc1 = "q48")
name_match_df <- tibble(rate = name_match_vec,
                        rate_us = names(name_match_vec))

## rename BT runs to match our names
fix_names <- function(x) {
        (x
            |> full_join(name_match_df, by = "rate")
            |> select(-rate)
            |> rename(rate = "rate_us")
        )
}

#### CONTRAST FUNCTIONS ####
## computing contrasts - from rates to geometric mean ratios
## rescued from 4b6dce5 on old_bayestest branch (!!)
gmean <- function(x, y) sqrt(x*y)  ## equivalent to exp((log(x) + log(y))/2)
cfun_nonlog <- function(q1, q2, q3, q4) {
    gmean(q1,q2)/gmean(q3,q4)
    ## or exp(
    ##        (log(q1) + log(q2))/2 -
    ##        (log(q3) + log(q4))/2
    ##    )
}

## function for getting contrasts from rates
get_contrasts <- function(results) {
  results$Log$results %>%
  mutate(gain_care_effect =   cfun_nonlog(q37,q48,q15,q26),
         gain_spawn_effect =  cfun_nonlog(q26,q48,q15,q37),
         gain_interaction =   cfun_nonlog(q15,q48,q37,q26),
         loss_care_effect =   cfun_nonlog(q73,q84,q51,q62),
         loss_spawn_effect =  cfun_nonlog(q62,q84,q51,q73),
         loss_interaction =   cfun_nonlog(q51,q84,q73,q62))
}

## read in three different formats (raw rates; rates as list of chains; contrasts)
bdir <- "bayestraits"
read_chains <- function(fn)      { readRDS(file.path(bdir, fn)) |> get_chains() }
read_chains_list <- function(fn) { readRDS(file.path(bdir, fn)) |> get_chains("mcmc.list") }
read_contrasts <- function(fn)   { readRDS(file.path(bdir, fn)) |> get_contrasts() }

if (!interactive()) pdf("bayestraits-pix.pdf", width = 16, height = 10)

## FIXME/TODO:
##  Compare with our results/side-by-side plots


prior_fac <- 1/1212  ## derived from scaling ratio of our vs BayesTraits trees (see bayestraits-run.R)

all_bt_model_results <- list.files(path=bdir,
                                pattern = "bt_model_(no)?data_(reg|rj)_(default|prior)")
names(all_bt_model_results) <- gsub("(bt_model_|\\.rds)", "", all_bt_model_results)
## should have a 2 x 2 x 2 result ({data, nodata} * {reg, rj} * {default, priors}) + 2
## we know 'short' RJ results are no good, drop them
all_bt_model_results <- all_bt_model_results[grep("rj_(priors|default)", names(all_bt_model_results), invert = TRUE)]
stopifnot(length(all_bt_model_results) == 6L)

no_rj <- grep("_reg_", names(all_bt_model_results))

## full info: parameters + contrasts
all_bt_results <- purrr::map_dfr(all_bt_model_results[no_rj], read_contrasts, .id = "model_run")
## q* parameters only
all_bt_chains <- purrr::map_dfr(all_bt_model_results, read_chains, .id = "model_run")
## diagnostics for each model

all_bt_diag <- (purrr::map_dfr(all_bt_model_results,
                           function(x) { read_chains_list(x) |> diagnostic_posterior() },
                           .id = "model_run")
    |> rename(rate = "Parameter")
    |> fix_names()
)

### load previous results

## what do we need to do to match BT and our results?
##  (1) match names of contrasts/rates
##  (2) rescale rates/contrasts

tar_load(ag_priorsamp)
tar_load(ag_mcmc_tb)
tar_load(contr_long_ag_mcmc_tb)
tar_load(contr_long_ag_priorsamp)


pivot_ours <- function(mcmc, nm = "ours") {
    (as.mcmc(mcmc) 
        |> as_tibble()
        |> mutate(iteration = seq(n()))
        |> pivot_longer(-iteration, names_to = "rate")
    )
}
add_info <- function(x,data = "data") {
    x |>
        mutate(data = data, method = "ours", priors = "priors", .before = 1)
}

our_chains <-    (map_dfr(ag_mcmc_tb, pivot_ours, .id = "chain")
    |> add_info()
    |> mutate(across(chain, as.numeric))
)

our_priors <-    (map_dfr(ag_priorsamp, pivot_ours, .id = "chain")
    |> add_info(data = "nodata")
    |> mutate(across(chain, as.numeric))
)

proc_contrasts <- function(x, data = "data") {
    (x
        |> add_info(data = data)
        |> filter(contrast != "intercept")
        |> filter(!grepl("netgain", rate))
        |> mutate(across(contrast,
                         ~ case_when(. == "pc" ~ "care_effect",
                                     . == "sc" ~ "spawn_effect",
                                     . == "pcxsc" ~ "interaction")),
                  name = sprintf("%s_%s", rate, contrast))
        |> mutate(across(value, exp))
    )
}

our_cpriors <- proc_contrasts(contr_long_ag_priorsamp, data = "nodata")
our_contrasts <- proc_contrasts(contr_long_ag_mcmc_tb)

## plot diagnostics
all_bt_diag_L <- (
    all_bt_diag
    |> pivot_longer(-c(model_run, rate))
    |> separate(model_run, into = c("data", "method", "priors"), sep = "_")
    |> mutate(across(priors, \(x) stringr::str_replace(x, "prior2-long", "priors")))
)

diag_plot <- all_bt_diag_L |> ggplot(aes(x = value, y = rate, colour = method, shape = interaction(priors, data))) +
    geom_point(size=5) +
    facet_wrap(~name, scale = "free") +
    scale_shape_manual(values = c(1, 16, 2, 17))  + ## open/closed x round/triangle
    scale_colour_manual(values = c("black", "red"))

print(diag_plot + labs(title = "diagnostics for all runs (reg + RJ long)"))

## results basically make sense
## ESS for RJ runs are high because they're **very long**. R-hat for q26 is bad
##   because it's practically always zero (this is ag0-1_pc0_sc1, gain of accessory
##  glands when 'no male parental care' + 'group spawning' (maybe very rare state?)

## get rid of RJ so we can focus on non-RJ
if (FALSE) {
    print(diag_plot %+% filter(all_bt_diag_L, method == "reg") +
          labs(title = "diagnostics for 'reg' runs only"))
}
      
## * nodata [triangles] gives high MCSE across the board (higher for default priors than ours)
## * worst R-hats are for priors.nodata (we probably don't care)

## trace plots for raw rate parameters ('q')
chains_long <- (all_bt_chains
    |> pivot_longer(starts_with("q"), names_to = "rate")
    |> fix_names()
)

gg_chains_1 <- ggplot(chains_long, aes(Iteration, value+min_val, colour = factor(chain))) +
    facet_grid(rate ~ model_run, scale = "free") + geom_line() +
    scale_y_log10()

print(gg_chains_1 + labs(title = "trace plots for raw rates ('q' params)") +
      ## https://stackoverflow.com/questions/48892826/rotate-strip-text-in-ggplot2
      theme(strip.text.y.right = element_text(angle = 0)))

## 
chains_long_2 <- (chains_long
    |> separate(model_run, into = c("data", "method", "priors"), sep = "_", remove = FALSE)
    |> mutate(across(priors, \(x) stringr::str_replace(x, "prior2-long", "priors")))
)

## log scale, but add min-value; narrower bandwidth for violin plots (adjust = 0.3)
gg_rate_violins <- ggplot(chains_long_2, aes(x = min_val + value, y = interaction(priors, data))) +
    geom_violin(aes(fill = factor(data)), alpha = 0.5, adjust = 0.4) +
    facet_grid(method~rate, scale = "free_y", space = "free") +
    scale_x_log10() +
    theme(panel.spacing = grid::unit(0, "lines")) +
    scale_fill_manual(values = c("green", "blue")) +
    labs(x="rate (+ min value)",
         y = "") +
    geom_vline(xintercept = 1.0, lty = 2)

chains_long_3 <- (
    chains_long_2
    ## |> mutate(across(value, ~ . / prior_fac))
    |> bind_rows(our_chains)
    |> bind_rows(our_priors)
)

sum_chains <- (chains_long_3
    |> filter(data == "data")
    |> group_by(method, priors, rate)
    |> summarise(across(value, mean))
    |> arrange(rate)
    |> group_by(rate)
    |> mutate(across(value, ~ .[method == "ours"]/.))
    |> filter(method == "reg")
)
print(sum_chains)
## BT mean rates are 2-20 times ours, not 1/prior_fac (1212x) ?
## maybe we can just compare p-values??

## FIXME: reverse legend order?
## divide into two blocks (gain/loss) for legibility
title_head <- "prior and posterior distributions for all rates & runs:"
print(gg_rate_violins %+% filter(chains_long_3, grepl("gain", rate)) +
      labs(title = paste(title_head, "GAINS")))
print(gg_rate_violins %+% filter(chains_long_3, grepl("loss", rate)) +
      labs(title = paste(title_head, "LOSSES")))

## contrasts: only for non-RJ examples
pivot_contrasts <- function(x) {
    (x
        |> select(c(any_of(c("model_run", "chain", "Iteration")),
                    matches("(interaction|effect)$")))
        |> pivot_longer(matches("(interaction|effect)"))
    )
}

chains_long_c <- (all_bt_results
    |> pivot_contrasts()
    |> separate(model_run, into = c("data", "method", "priors"), remove = FALSE)
)

## trace plots of contrasts
print(gg_chains_1 %+% (chains_long_c |> rename(rate = "name")) + labs(title = "trace plots for contrasts (log scale)"))
## without log scale
## print(chains_1 %+% chains_long_c + scale_y_continuous() + labs(title = "trace plots for contrasts (non-log scale)"))

chains_long_c_2 <- (
    bind_rows(chains_long_c, our_contrasts, our_cpriors)
    |> mutate(across(data, ~factor(., levels = c("nodata", "data"))))
)
## show all posterior distributions
ggplot(chains_long_c_2, aes(x = value, y = interaction(priors, data))) +
    geom_violin(aes(fill = factor(data)), alpha = 0.5) +
    facet_grid(method~name, scale = "free_y", space = "free") +
    scale_x_log10() +
    theme(panel.spacing = grid::unit(0, "lines")) +
    scale_fill_manual(values = c("green", "blue")) +
    labs(title = "posterior distributions for all contrasts & runs") +
    geom_vline(xintercept = 1.0, lty = 2)
                              
## confidence intervals and p-values

Signif <- function(pv) symnum(pv, corr = FALSE, na = FALSE, 
                 cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 symbols = c("***", "**", "*", ".", " "))

csum <- (chains_long_c_2
    |> group_by(name,method, data)
    |> summarise(
           lwr = quantile(value, 0.025),
           lwr50 = quantile(value, 0.25),
           med = median(value),
           upr = quantile(value, 0.975),
           upr50 = quantile(value, 0.75),
           pval = 2*min(mean(value<1), mean(value>1)),
           psym = as.character(Signif(pval)),
           .groups = "drop")
)

(csum
    |> select(name, method, data, pval)
    |> pivot_wider(names_from = method, values_from = pval)
)


pd <- position_dodge(width = 0.25)
ggplot(csum, aes(y = name, x = med, colour = method)) +
    geom_pointrange(aes(xmin = lwr, xmax = upr),
                    position = pd) +
    geom_pointrange(aes(xmin = lwr50, xmax = upr50), linewidth = 2,
                    position = pd) +
    scale_x_log10() +
    facet_wrap(~data) +
    scale_colour_manual(values = c("darkgreen", "blue")) +
    geom_vline(xintercept = 1, lty = 2) +
    expand_limits(x=100) +
    geom_text(aes(x = 100, label = psym), size = 10, position = position_dodge(width = 0.5),
              show.legend = FALSE) +
    labs(x="contrast", y = "", title = "50% and 95% CIs for contrasts (prior and posterior)")

## raftery.diag(chains1)  ## suggests we have to run longer (~4000 samples == 8 x 500)
## g1 <- geweke.diag(chains1)
## calculate 2-tailed p-values for geweke diagnostics based on the Z scores
## sort(2*pnorm(abs(g1$z), lower.tail = FALSE))

##### WEIGHTING ####

## we simulate character histories with simmap, supplying our model
## from that we get the amount of time spent in each state (1-8)
## we call the time spent in state 1 "time1" and so on

## g mean written the long way:
## exp( ((log(q12)*(1/2)) + (log(q34)*(1/2))) )
## sub first (1/2) to (time1/(time1+time3)) 

## function to get weighted mean
gmean_weighted <- function(q12, q34, time1, time3) {
  exp(
    (log(q12)*
    (time1/(time1+time3))) +
    (log(q34)*
    (time3/(time1+time3))) 
     )
}

## weighted contrasts
cfun_nonlog_weighted <- function(q12, q34, q56, q78, time1, time3, time5, time7) {
    gmean(q12,q34,time1,time3)/gmean(q56,q78,time5,time7)
}

## function for getting weighted contrasts from rates
get_contrasts_weighted <- function(results) {
  results$Log$results %>% 
    mutate(gain_care_effect =   cfun_nonlog(q37,q48,q15,q26,time3,time4,time1,time2),
           gain_spawn_effect =  cfun_nonlog(q26,q48,q15,q37,time2,time4,time1,time3),
           gain_interaction =   cfun_nonlog(q15,q48,q37,q26,time1,time4,time3,time2),
           loss_care_effect =   cfun_nonlog(q73,q84,q51,q62,time7,time8,time5,time6),
           loss_spawn_effect =  cfun_nonlog(q62,q84,q51,q73,time6,time8,time5,time7),
           loss_interaction =   cfun_nonlog(q51,q84,q73,q62,time5,time8,time7,time6))
}

if (!interactive()) dev.off()


## what do we want? prior_contr_ci, ag_contr_gainloss
q()


#### RATE DESCRIPTIONS ####
## q12 = spawnGain
## q13 = careGain
## q14 = 0
## q15 = ag0-1_pc0_sc0
## q16 = 0
## q17 = 0
## q18 = 0
## q21 = spawnLoss
## q23 = 0
## q24 = careGain
## q25 = 0
## q26 = ag0-1_pc0_sc1
## q27 = 0
## q28 = 0
## q31 = careLoss
## q32 = 0
## q34 = spawnGain
## q35 = 0
## q36 = 0
## q37 = ag0-1_pc1_sc0
## q38 = 0
## q41 = 0 
## q42 = careLoss
## q43 = spawnLoss
## q45 = 0
## q46 = 0
## q47 = 0
## q48 = ag0-1_pc1_sc1
## q51 = ag1-0_pc0_sc0
## q52 = 0
## q53 = 0
## q54 = 0
## q56 = spawnGain
## q57 = careGain
## q58 = 0
## q61 = 0
## q62 = ag1-0_pc0_sc1
## q63 = 0
## q64 = 0
## q65 = spawnLoss
## q67 = 0
## q68 = careGain
## q71 = 0
## q72 = 0
## q73 = ag1-0_pc1_sc0
## q74 = 0
## q75 = careLoss
## q76 = 0
## q78 = spawnGain
## q81 = 0
## q82 = 0
## q83 = 0
## q84 = ag1-0_pc1_sc1
## q85 = 0

## are these duplicates/restricted??
## q86 = careLoss
## q87 = spawnLoss

