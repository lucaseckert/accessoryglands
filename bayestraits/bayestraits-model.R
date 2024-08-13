#### BayesTraits Model analysis ####

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

## computing contrasts
## FIXME: we need to be more careful about computing contrasts,
##  i.e. on the *log scale* we were computing something (q1 + q2)/2 - (q3 + q4)/2 = m1 - m2
##  when we exponentiate we get exp(m1)/exp(m2) -- so the ratio is correct
##  but we should be taking the _geometric mean_ of the rates, not the arithmetic mean
##  e.g.
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

#### DATA, REGULAR, DEFAULT PRIORS ####

read_contrasts <- function(fn) { readRDS(file.path("bayestraits", fn)) |> get_contrasts() }
read_chains <- function(fn) { readRDS(file.path("bayestraits", fn)) |> get_chains() }
read_chains_list <- function(fn) { readRDS(file.path("bayestraits", fn)) |> get_chains("mcmc.list") }

all_model_results <- list.files(path="bayestraits", pattern = "bt_model_(no)?data_.*")
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

chains_long_c <- (all_results
    |> select(c("model_run", "chain", "Iteration", matches("(interaction|effect)$")))
    |> pivot_longer(matches("(interaction|effect)"))
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
q()


#### OLD CODE ####

# cfun_log <- function(q1, q2, q3, q4) {
#     (q1 + q2)/2 - (q3 + q4)/2
# }

## LUCAS: decide how to fix this. 
##   e.g.

## FIX ME !!!

## either (1) convert rates to log-rates, i.e.
## rates[,rate_cols] <- log(rates[,rate_cols])

##and use code below
##  (or change everything to cfun_log for compactness)
## OR (2) use cfun_nonlog below

# get_contrasts <- function(rates) {
#     mutate(rates,
#            gain_care_effect =   (q37+q48)/2-(q15+q26)/2,
#            gain_spawn_effect =  (q26+q48)/2-(q15+q37)/2,
#            gain_interaction =   (q15+q48)/2-(q37+q26)/2,
#            loss_care_effect =   (q73+q84)/2-(q51+q62)/2,
#            loss_spawn_effect =  (q62+q84)/2-(q51+q73)/2,
#            loss_interaction =   (q51+q84)/2-(q73+q62)/2)
# }

## gain contrasts
gg_gain <- dplyr::select(contrasts, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin(fill="gray")+
  theme_bw()+
  scale_x_continuous(trans = "log10", limits = c(1e-2, 1e2))

print(gg_gain) + labs(title = "gain contrasts, our priors, non-RJ")
## main issue here: gain_spawn_effect positive rather than negative?

## loss contrasts
cdat <- dplyr::select(contrasts, loss_care_effect:loss_interaction) %>% 
  pivot_longer(cols=loss_care_effect:loss_interaction, names_to = "contrast")

gg_loss <- gg_gain %+% cdat
print(gg_loss)

## checking out the acceptance rate
summary(schedule)

## command vector
command_vec_rj <- get_command(prior = "RevJump uniform 0 100")
#results_rj<-bayestraits(data_dataless, trees, command_vec_rj)
#saveRDS(results_rj, file = "bayestraits/bt_model_rj.rds")

## reading in results
results_rj<-readRDS("bayestraits/bt_model_rj.rds")
rates_rj<-results_rj$Log$results
summary(rates_rj)

## computing contrasts
contrasts_rj <- get_contrasts(rates_rj)
cdef_rj <- dplyr::select(contrasts_rj, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast")

gg_gain %+% cdef_rj

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
## q86 = careLoss
## q87 = spawnLoss



