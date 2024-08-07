#### BayesTraits Model ####

## run from head directory of repo
## (for target loading, readRDS)

## packages
library(btw)
library(targets)
library(tidyverse)
library(ape)
library(coda)
library(bayestestR)
## see http://bbolker.github.io/bbmisc/bayes/examples.html
## for convergence diagnostics etc.
## however, we would have to make BT output look like results from
##  a stan fit/write appropriate methods ...

if (!interactive()) pdf("bayestraits-pix.pdf")
## loading trees
tar_load(treeblock)
trees <- do.call(c, treeblock)
trees <- .compressTipLabel(trees)

## loading data

## FIXME/TODO: unify ggplots
##  Compare with our results/side-by-side plots
##  Bayes diagnostics (R-hat, trace plots ...)


## FIXME: can we automate this?
## (4*ag + 2*pc + sc will work for the non-missing cases ...)
tar_load(ag_compdata_tb)

data <- ag_compdata_tb$data %>%
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
    select(species,state)

summary(data)

#### PRIOR SCALING ####

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


#### COMMAND FUNCTION ####

##  q14 q16 q17 q18 q23 q25 q27 q28 q32 q35 q36 q38 q41 q45 q46 q47 q52 q53 q54 q58 q61 q63 q64 q67 q71 q72 q74 q76 q81 q82 q83 q85 0", ## impossible rates 
zero_rates <- c(14, 16:18, 23, 25, 27:28, 32, 35:36, 38, 41, 45:47, 52:54,
                58, 61, 63:64, 67, 71:72, 74, 76, 81:83, 85)
qz <- paste("q", zero_rates, sep ="", collapse = " ")
bt_command <- function(prior = NULL, iterations = 51e4, burnin = 1e4, seed = 101) {
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
               sprintf("Seed %d", seed)
               )
    return(cvec)
}


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
plot_contrasts<-function(contrasts) {
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


#### DATA, REGULAR, DEFAULT PRIOS ####

## command vector
command_vec_data_reg_default<- bt_command(prior = NULL)
#results_data_reg_default<-bayestraits(data, trees, command_vec_data_reg_default)
#saveRDS(results_data_reg_default, file = "bayestraits/bt_model_data_reg_default.rds")

## reading in results
results_data_reg_default <- readRDS("bayestraits/bt_model_data_reg_default.rds")

## computing contrasts
contrasts_data_reg_default<-get_contrasts(results_data_reg_default)

## plotting
plot_contrasts(contrasts_data_reg_default)


##### DATA, REGULAR, OUR PRIORS ####

## command vector
command_vec_data_reg_priors<- bt_command(prior = prior1)
#results_data_reg_priors<-bayestraits(data, trees, command_vec_data_reg_priors)
#saveRDS(results_data_reg_priors, file = "bayestraits/bt_model_data_reg_priors.rds")

## read in results
results_data_reg_priors<-readRDS("bayestraits/bt_model_data_reg_priors.rds")

## computing contrasts
contrasts_data_reg_priors<-get_contrasts(results_data_reg_priors)

## plotting
plot_contrasts(contrasts_data_reg_priors)


#### DATA-LESS, REGULAR, DEFAULT PRIORS ####

## data-less data
data_dataless <- mutate(data, state=12345678)
summary(data_dataless)

## command vector
command_vec_nodata_reg_default<-bt_command(prior = NULL)
#results_nodata_reg_default<-bayestraits(data_dataless, trees, command_vec_nodata_reg_default)
#saveRDS(results_nodata_reg_default, file = "bayestraits/bt_model_nodata_reg_default.rds")

## reading results
results_nodata_reg_default<-readRDS("bayestraits/bt_model_nodata_reg_default.rds")

## contrasts
contrasts_nodata_reg_default<-get_contrasts(results_nodata_reg_default)

## plotting
plot_contrasts(contrasts_nodata_reg_default)


#### DATA-LESS, REGULAR, OUR PRIORS ####

## command vector
command_vec_nodata_reg_priors<-bt_command(prior = prior1)
#results_nodata_reg_priors<-bayestraits(data_dataless, trees, command_vec_nodata_reg_priors)
#saveRDS(results_nodata_reg_priors, file = "bayestraits/bt_model_nodata_reg_priors.rds")

## reading in results
results_nodata_reg_priors<-readRDS("bayestraits/bt_model_nodata_reg_priors.rds")

## computing contrasts
contrasts_nodata_reg_priors<-get_contrasts(results_nodata_reg_priors)

## plotting
plot_contrasts(contrasts_nodata_reg_priors)


#### DATA, RJ, DEFAULT ####

## command vector
command_vec_data_rj_default<-bt_command(prior = "RevJump uniform 0 100")
#results_data_rj_default<-bayestraits(data, trees, command_vec_data_rj_default)
#saveRDS(results_data_rj_default, file = "bayestraits/bt_model_data_rj_default.rds")

## reading in results
results_data_rj_default<-readRDS("bayestraits/bt_model_data_rj_default.rds")

## computing contrasts
contrasts_data_rj_default<-get_contrasts(results_data_rj_default)

## plotting
plot_contrasts(contrasts_data_rj_default)

#### DATA-LESS, RJ, DEFAULT ####

## command vector
command_vec_nodata_rj_default<-bt_command(prior = "RevJump uniform 0 100")
#results_nodata_rj_default<-bayestraits(data_dataless, trees, command_vec_nodata_rj_default)
#saveRDS(results_nodata_rj_default, file = "bayestraits/bt_model_nodata_rj_default.rds")

## reading in results
results_nodata_rj_default<-readRDS("bayestraits/bt_model_nodata_rj_default.rds")

## computing contrasts
contrasts_nodata_rj_default<-get_contrasts(results_nodata_rj_default)

## plotting
plot_contrasts(contrasts_nodata_rj_default)


#### DATA, RJ, OUR PRIORS ####

## command vector
command_vec_data_rj_priors<-bt_command(prior = prior2)
#results_data_rj_priors<-bayestraits(data, trees, command_vec_data_rj_priors)
#saveRDS(results_data_rj_priors, file = "bayestraits/bt_model_data_rj_priors.rds")

## reading in results
results_data_rj_priors<-readRDS("bayestraits/bt_model_data_rj_priors.rds")

## computing contrasts
contrasts_data_rj_priors<-get_contrasts(results_data_rj_priors)

## plotting
plot_contrasts(contrasts_data_rj_priors)


#### DATA-LESS, RJ, OUR PRIORS ####

## command vector
command_vec_nodata_rj_priors<-bt_command(prior = prior2)
#results_nodata_rj_priors<-bayestraits(data_dataless, trees, command_vec_nodata_rj_priors)
#saveRDS(results_nodata_rj_priors, file = "bayestraits/bt_model_nodata_rj_priors.rds")

## reading in results
results_nodata_rj_priors<-readRDS("bayestraits/bt_model_nodata_rj_priors.rds")

## computing contrasts
contrasts_nodata_rj_priors<-get_contrasts(results_nodata_rj_priors)

## plotting
plot_contrasts(contrasts_nodata_rj_priors)


#### DIAGNOSTICS ####

## run in parallel???
options(bt_path= "BayesTraitsV4.0.0-Linux")
options(bt_bin = "BayesTraitsV4")

## FIXME: cache this more sensibly
## FIXME: does this work OK with BT V4?
## FIXME: split this file
## FIXME: run bayestraits multiple times with distinct seeds
## (and possibly shorter chains?) so that
##  we can calculate Gelman-Rubin diagnostics (R-hat)

## reading in results
results <- readRDS("bayestraits/bt_model_demo.rds")
rates <- results$Log$results
options <- results$Log$options
schedule <- results$Schedule$header

## 4:59
rate_cols <- grep("^q[0-9]+", colnames(results$Log$results))

get_chains <- function(results) {
  chains <- as.mcmc(results$Log$results[,rate_cols])  ## q** values only
  cols_disallowed <- which(apply(chains==0, 2, all)) ## forbidden/boring
  dupes <- c("q24","q57","q68", ## care gain
             "q42","q75","q86", ## care loss
             "q34","q56","q78", ## spawn gain
             "q43","q65","q87" ## spawn loss
  )
  cols_dupes <- match(dupes, colnames(chains))
  
  return(chains[, -c(cols_dupes, cols_disallowed)])
}

chains1 <- get_chains(results)
lattice::xyplot(chains1, aspect = "fill", layout = c(4,3),
                scales = list(y = list(log = 10)))


raftery.diag(chains1)  ## suggests we have to run longer (~4000 samples == 8 x 500)
g1 <- geweke.diag(chains1)
## calculate 2-tailed p-values for geweke diagnostics based on the Z scores
sort(2*pnorm(abs(g1$z), lower.tail = FALSE))


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
gg_gain <- select(contrasts, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin(fill="gray")+
  theme_bw()+
  scale_x_continuous(trans = "log10", limits = c(1e-2, 1e2))

print(gg_gain) + labs(title = "gain contrasts, our priors, non-RJ")
## main issue here: gain_spawn_effect positive rather than negative?

## loss contrasts
cdat <- select(contrasts, loss_care_effect:loss_interaction) %>% 
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
cdef_rj <- select(contrasts_rj, gain_care_effect:gain_interaction) %>% 
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


if (!interactive()) dev.off()
