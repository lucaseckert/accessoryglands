#### BayesTraits Model ####

## run from ???

## packages
library(btw)
library(targets)
library(tidyverse)
library(ape)

## loading trees
tar_load(treeblock)
trees <- do.call(c, treeblock)
trees <- .compressTipLabel(trees)

## loading data
tar_load(ag_compdata_tb)
data<-ag_compdata_tb$data %>% mutate(state=factor(case_when((ag==0 & pc==0 & sc==0) ~ 1,
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

#### REAL MODEL ####

## command vector
command_vec<- c("1", ## MultiState
                "2", ## MCMC
                "ScaleTrees", ## scaling branch lengths to a mean of 0.1
                "AddTag Root Erpetoichthys_calabaricus Mugil_liza", ## adding a tag at the root
                "Fossil Node01 Root 2", ## fossilizing the root
                "Res q14 q16 q17 q18 q23 q25 q27 q28 q32 q35 q36 q38 q41 q45 q46 q47 q52 q53 q54 q58 q61 q63 q64 q67 q71 q72 q74 q76 q81 q82 q83 q85 0", ## impossible rates 
                "Res q13 q24 q57 q68", ## care gain
                "Res q31 q42 q75 q86", ## care loss
                "Res q12 q34 q56 q78", ## spawn gain
                "Res q21 q43 q65 q87", ## spawn loss
                "PriorAll lognormal 4.236 1.41", ## setting priors lognormal with mean and sd
                "Iterations 510000", ## iterations, including burn-in
                "Burnin 10000") ## burn-in

## if you want to run it again 
<<<<<<< HEAD
#results<-bayestraits(data, trees, command_vec)
#saveRDS(results, file = "bayestraits/bt_model_demo.rds")
=======
## results <- bayestraits(data, trees, command_vec)
## saveRDS(results, file = "bayestraits/bt_model_demo.rds")
>>>>>>> ea75bbbcb6502d8cbe26a38457fe1ffa1e27e1f5

## reading in results
results <- readRDS("bayestraits/bt_model_demo.rds")
rates <- results$Log$results
options <- results$Log$options
schedule <- results$Schedule$header

## computing contrasts
contrasts<-mutate(rates, gain_care_effect = ((q37+q48)/2)/((q15+q26)/2),
                         gain_spawn_effect = ((q26+q48)/2)/((q15+q37)/2),
                         gain_interaction = ((q15+q48)/2)/((q37+q26)/2),
                         loss_care_effect = ((q73+q84)/2)/((q51+q62)/2),
                         loss_spawn_effect = ((q62+q84)/2)/((q51+q73)/2),
                         loss_interaction = ((q51+q84)/2)/((q73+q62)/2))

## gain contrasts
select(contrasts, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin()+
  theme_bw()+
  scale_x_continuous(trans = "log10")

## loss contrasts
select(contrasts, loss_care_effect:loss_interaction) %>% 
  pivot_longer(cols=loss_care_effect:loss_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin(fill = "gray")+
  theme_bw()+
  scale_x_continuous(trans = "log10")

## checking out the acceptance rate
summary(schedule)

<<<<<<< HEAD
#### MODEL W/ DEFAULT PRIORS ####

## command vector
command_vec_default<- c("1", ## MultiState
                         "2", ## MCMC
                         "ScaleTrees", ## scaling branch lengths to a mean of 0.1
                         "AddTag Root Erpetoichthys_calabaricus Mugil_liza", ## adding a tag at the root
                         "Fossil Node01 Root 2", ## fossilizing the root
                         "Res q14 q16 q17 q18 q23 q25 q27 q28 q32 q35 q36 q38 q41 q45 q46 q47 q52 q53 q54 q58 q61 q63 q64 q67 q71 q72 q74 q76 q81 q82 q83 q85 0", ## impossible rates 
                         "Res q13 q24 q57 q68", ## care gain
                         "Res q31 q42 q75 q86", ## care loss
                         "Res q12 q34 q56 q78", ## spawn gain
                         "Res q21 q43 q65 q87", ## spawn loss
                         "Iterations 510000", ## iterations, including burn-in
                         "Burnin 10000") ## burn-in
#results_default<-bayestraits(data, trees, command_vec_default)
#saveRDS(results_default, file = "bayestraits/bt_model_default.rds")

## reading in results
results_default<-readRDS("bayestraits/bt_model_default.rds")
rates_default<-results_default$Log$results
summary(rates_default)

## computing contrasts
contrasts_default<-mutate(rates_default, gain_care_effect = ((q37+q48)/2)/((q15+q26)/2),
                  gain_spawn_effect = ((q26+q48)/2)/((q15+q37)/2),
                  gain_interaction = ((q15+q48)/2)/((q37+q26)/2),
                  loss_care_effect = ((q73+q84)/2)/((q51+q62)/2),
                  loss_spawn_effect = ((q62+q84)/2)/((q51+q73)/2),
                  loss_interaction = ((q51+q84)/2)/((q73+q62)/2))
select(contrasts_default, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin()+
  theme_bw()+
  scale_x_continuous(trans = "log10")

#### DATA-LESS MODEL W/ OUR PRIORS ####

data_dataless<-mutate(data, state=12345678)
summary(data_dataless)

## command vector
command_vec_dataless<- c("1", ## MultiState
                         "2", ## MCMC
                         "ScaleTrees", ## scaling branch lengths to a mean of 0.1
                         "AddTag Root Erpetoichthys_calabaricus Mugil_liza", ## adding a tag at the root
                         "Fossil Node01 Root 2", ## fossilizing the root
                         "Res q14 q16 q17 q18 q23 q25 q27 q28 q32 q35 q36 q38 q41 q45 q46 q47 q52 q53 q54 q58 q61 q63 q64 q67 q71 q72 q74 q76 q81 q82 q83 q85 0", ## impossible rates 
                         "Res q13 q24 q57 q68", ## care gain
                         "Res q31 q42 q75 q86", ## care loss
                         "Res q12 q34 q56 q78", ## spawn gain
                         "Res q21 q43 q65 q87", ## spawn loss
                         "PriorAll lognormal 4.236 1.41", ## setting priors lognormal with mean and sd
                         "Iterations 510000", ## iterations, including burn-in
                         "Burnin 10000") ## burn-in
#results_dataless<-bayestraits(data_dataless, trees, command_vec_dataless)
#saveRDS(results_dataless, file = "bayestraits/bt_model_dataless.rds")

## reading in results
results_dataless<-readRDS("bayestraits/bt_model_dataless.rds")
rates_dataless<-results_dataless$Log$results
summary(rates_dataless)

## computing contrasts
contrasts_dataless<-mutate(rates_dataless, gain_care_effect = ((q37+q48)/2)/((q15+q26)/2),
                          gain_spawn_effect = ((q26+q48)/2)/((q15+q37)/2),
                          gain_interaction = ((q15+q48)/2)/((q37+q26)/2),
                          loss_care_effect = ((q73+q84)/2)/((q51+q62)/2),
                          loss_spawn_effect = ((q62+q84)/2)/((q51+q73)/2),
                          loss_interaction = ((q51+q84)/2)/((q73+q62)/2))
select(contrasts_dataless, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin()+
  theme_bw()+
  scale_x_continuous(trans = "log10")

#### RJ DATA-LESS MODEL W/ DEFAULT PRIORS ####

## command vector
command_vec_rj<- c("1", ## MultiState
                   "2", ## MCMC
                   "ScaleTrees", ## scaling branch lengths to a mean of 0.1
                   "AddTag Root Erpetoichthys_calabaricus Mugil_liza", ## adding a tag at the root
                   "Fossil Node01 Root 2", ## fossilizing the root
                   "Res q14 q16 q17 q18 q23 q25 q27 q28 q32 q35 q36 q38 q41 q45 q46 q47 q52 q53 q54 q58 q61 q63 q64 q67 q71 q72 q74 q76 q81 q82 q83 q85 0", ## impossible rates 
                   "Res q13 q24 q57 q68", ## care gain
                   "Res q31 q42 q75 q86", ## care loss
                   "Res q12 q34 q56 q78", ## spawn gain
                   "Res q21 q43 q65 q87", ## spawn loss
                   "RevJump uniform 0 100", ## setting priors 
                   "Iterations 510000", ## iterations, including burn-in
                   "Burnin 10000") ## burn-in
#results_rj<-bayestraits(data_dataless, trees, command_vec_rj)
#saveRDS(results_rj, file = "bayestraits/bt_model_rj.rds")

## reading in results
results_rj<-readRDS("bayestraits/bt_model_rj.rds")
rates_rj<-results_rj$Log$results
summary(rates_rj)

## computing contrasts
contrasts_rj<-mutate(rates_rj, gain_care_effect = ((q37+q48)/2)/((q15+q26)/2),
                          gain_spawn_effect = ((q26+q48)/2)/((q15+q37)/2),
                          gain_interaction = ((q15+q48)/2)/((q37+q26)/2),
                          loss_care_effect = ((q73+q84)/2)/((q51+q62)/2),
                          loss_spawn_effect = ((q62+q84)/2)/((q51+q73)/2),
                          loss_interaction = ((q51+q84)/2)/((q73+q62)/2))
select(contrasts_rj, gain_care_effect:gain_interaction) %>% 
  pivot_longer(cols=gain_care_effect:gain_interaction, names_to = "contrast") %>% 
  ggplot(aes(x=value, y=contrast))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_violin()+
  theme_bw()+
  scale_x_continuous(trans = "log10")

=======
>>>>>>> ea75bbbcb6502d8cbe26a38457fe1ffa1e27e1f5
#### Rate Descriptions ####
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


