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

source("R/utils.R")

if (!interactive()) pdf("bayestraits-pix.pdf", width = 16, height = 10)

## FIXME/TODO: unify ggplots
##  Compare with our results/side-by-side plots
##  Bayes diagnostics (R-hat, trace plots ...)

#### CONTRAST FUNCTIONS ####

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
get_contrasts <- function(results) {
  results$Log$results %>% 
  mutate(gain_care_effect =   cfun_nonlog(q37,q48,q15,q26),
         gain_spawn_effect =  cfun_nonlog(q26,q48,q15,q37),
         gain_interaction =   cfun_nonlog(q15,q48,q37,q26),
         loss_care_effect =   cfun_nonlog(q73,q84,q51,q62),
         loss_spawn_effect =  cfun_nonlog(q62,q84,q51,q73),
         loss_interaction =   cfun_nonlog(q51,q84,q73,q62))
}

## function for plotting contrasts
plot_contrasts <- function(contrasts) {
  df_name <- deparse(substitute(contrasts))
  contrasts %>% 
  pivot_longer(cols=gain_care_effect:loss_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin(fill="gray")+
  scale_x_continuous(trans = "log10")+
  theme_bw()+
  ggtitle(df_name)
}


prior_fac <- 1/1212  ## derived 
#### DATA, REGULAR, DEFAULT PRIORS ####

read_contrasts <- function(fn) { readRDS(file.path("bayestraits", fn)) |> get_contrasts() }
read_chains <- function(fn) { readRDS(file.path("bayestraits", fn)) |> get_chains() }
read_chains_list <- function(fn) { readRDS(file.path("bayestraits", fn)) |> get_chains("mcmc.list") }

all_model_results <- list.files(path="bayestraits", pattern = "bt_model_(no)?data_(reg|rj)_(default|priors).*")
## should have a 2 x 2 x 2 result ({data, nodata} * {reg, rj} * {default, priors})
stopifnot(length(all_model_results) == 8)

names(all_model_results) <- gsub("(bt_model_|\\.rds)", "", all_model_results)
## full info: parameters + contrasts
all_results <- purrr::map_dfr(all_model_results, read_contrasts, .id = "model_run")
## q* parameters only
all_chains <- purrr::map_dfr(all_model_results, read_chains, .id = "model_run")
## diagnostics for each model
all_diag <- purrr::map_dfr(all_model_results,
                           function(x) { read_chains_list(x) |> diagnostic_posterior() },
                           .id = "model_run")

## plot diagnostics
all_diag_L <- all_diag |> pivot_longer(-c(model_run, Parameter)) |>
    separate(model_run, into = c("data", "method", "priors"))

diag_plot <- all_diag_L |> ggplot(aes(x = value, y = Parameter, colour = method, shape = interaction(priors, data))) +
    geom_point(size=5) +
    facet_wrap(~name, scale = "free") +
    scale_shape_manual(values = c(1, 16, 2, 17))  + ## open/closed x round/triangle
    scale_colour_manual(values = c("black", "red"))

print(diag_plot + labs(title = "diagnostics for all runs (reg + RJ)"))

## results basically make sense
## * rj [red] is awful (low ESS, high R-hat most of the time) [trace plots will illustrate this even more strongly]

## get rid of RJ so we can focus on non-RJ
print(diag_plot %+% filter(all_diag_L, method == "reg") +
      labs(title = "diagnostics for 'reg' runs only"))
      

## * nodata [triangles] gives high MCSE across the board (higher for default priors than ours)
## * worst R-hats are for priors.nodata (we probably don't care)


## trace plots for raw rate parameters ('q')
chains_long <- (all_chains
    |> pivot_longer(starts_with("q"))
)
chains_1 <- ggplot(chains_long, aes(Iteration, value, colour = factor(chain))) +
    facet_grid(name ~ model_run, scale = "free_y") + geom_line() +
    scale_y_log10()
print(chains_1 + labs(title = "trace plots for raw rates ('q' params)"))

    |> pivot_longer(starts_with("q"))
)

pivot_contrasts <- function(x) {
    (x
        |> select(c(any_of(c("model_run", "chain", "Iteration")),
                    matches("(interaction|effect)$")))
        |> pivot_longer(matches("(interaction|effect)"))
    )
}

chains_long_c <- (all_results
    |> pivot_contrasts()
    |> separate(model_run, into = c("data", "method", "priors"), remove = FALSE)
)

## trace plots of contrasts
chains_1 %+% chains_long_c + labs(title = "trace plots for contrasts (log scale)")
## without log scale
chains_1 %+% chains_long_c + scale_y_continuous() + labs(title = "trace plots for contrasts (non-log scale)")

## show all posterior distributions

ggplot(chains_long_c, aes(x = value, y = interaction(priors, data))) + geom_violin(aes(fill = factor(data)), alpha = 0.5) +
    facet_grid(method~name) + scale_x_log10() + theme(panel.spacing = grid::unit(0, "lines")) +
    scale_fill_manual(values = c("green", "blue")) +
    labs(title = "posterior distributions for all contrasts & runs") +
    geom_vline(xintercept = 1.0, lty = 2)


## exploring experiments

## ex_fn <- "bt_model_data_rj_prior-exp10-long.rds"
ex_fn <- "bt_model_data_rj_prior2-long.rds"
ex_chain <- read_chains(ex_fn)
## thin a bit more
nn <- nrow(ex_chain)
ss <- seq(1, nrow(ex_chain), by = nn %/% 4800)
ex_chain <- ex_chain |> slice(ss)
ex_contrasts <- read_contrasts(ex_fn) |> slice(ss)

ex_chain_L <- (ex_chain
    ## |> slice(seq(1, nrow(ex_chain), by = 10))
    |> pivot_longer(starts_with("q"))
)

## mostly, sort of, OK except for q26
gg_excontrasts <- ggplot(ex_chain_L, aes(Iteration, value+1e-4)) +
    geom_line(aes(colour = factor(chain))) +
    scale_y_log10() +
    facet_wrap(~name, scale = "free") +
    labs(title = "trace plots for raw rates ('q' params)")
print(gg_excontrasts)

## just one ...
gg_excontrasts %+% filter(ex_chain_L, name == "q48")

gg_violin <- ggplot(ex_chain_L, aes(value+1e-4, name)) +
    geom_violin(aes(fill = factor(chain)), alpha = 0.5, position = "identity",
                bw = 0.05)
print(gg_violin)

## not quite what I wanted ...
gg_violin %+% filter(ex_chain_L, name == "q48")

filter(ex_chain_L, name == "q48") |>
    ggplot(aes(x=value)) +
    geom_histogram(bins = 100) +
    ## log10 doesn't really make sense but it shows what I want ...
    scale_y_log10()

filter(ex_chain_L, name == "q48") |>
    ggplot(aes(x=value)) +
    geom_density(bw=0.01) + 
    ## log10 doesn't really make sense but it shows what I want ...
    scale_y_log10(limits = c(1e-3, NA))

## geom-mean contrasts don't work well when rates are sometimes
##  set to zero
ggplot(pivot_contrasts(ex_contrasts),
       aes(value, name)) + geom_violin(fill = "gray")

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

##
tar_load(ag_mcmc_tb)
colnames(ag_mcmc_tb[[1]])
name_match_vec <- c(loss.sc="q21",
                loss.pc="q42",  ## also q86
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
name_match_df <- tibble(rate_bt = name_match_vec,
                        rate_us = names(name_match_vec))

## compare rates with RJ rates

tar_load(contr_long_ag_mcmc_tb)
tar_load(ag_mcmc_tb)
aa <- ag_mcmc_tb[[1]]
as.tibble.mcmc <- function(x) {
    x2 <- unclass(x)
    attr(x2, "mcpar") <- NULL
    x2 |> as_tibble(rownames = "iteration")
}

(                                   
    ag_mcmc_tb 
    |> purrr::map_dfr(as.tibble.mcmc, .id = "chain")
    |> tidyr::pivot_longer(-c(chain, iteration), names_to = "rate")
    |> full_join(name_match
)

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



