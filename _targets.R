library("targets")
library("tarchetypes")
library("tidyverse")
source("R/utils.R")
source("R/functions.R")
source("R/mcmc.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = pkgList)
redo_mcmc <- TRUE
run_mcmc <- function() {
  tar_cue(if (redo_mcmc) "thorough" else "never")
}

## TO DEBUG: set tar_option_set(debug = "target_name"); tar_make(callr_function = NULL)
## tar_option_set(debug = "ag_pcsc_pars")
## tar_option_set(debug = "ag_compdata")

## to play around interactively:
## library(targets); for (f in c("utils", "functions", "mcmc")) source(sprintf("R/%s.R", f))
## tar_load(everything())

## tar_plan is from tarchetypes::tar_plan()
data_input_targets <- tar_plan(
    ## trait data file
    tar_target(
        ag_binary_trait_file,
        "data/binaryTraitData.csv",
        format = "file"
    ),

    ## imputed phylogenies file
    tar_target(
        treeblock_file,
        "data/treeBlock.rds",
        format = "file"
    ),

    ## read trait data
    tar_target(
        full_ag_data,
        read_csv(ag_binary_trait_file, col_types = cols())
    ),

    ## read imputed phylogenies
    tar_target(
        treeblock,
        readRDS("data/treeBlock.rds") %>% purrr::map(scale_phylo)
    ),

    ## read full (genetic-data-only) phylogenies
    tar_target(
        fishtree_phylo,
        (suppressWarnings(
            fishtree_phylogeny(full_ag_data$species))
          %>% scale_phylo()
        )
    ),

    ## FIXME: use target_map() to combine these two
    ## read trait file, combine with fishtree phylo / trim
    tar_target(
        ag_compdata,
        get_ag_data(full_ag_data, phylo = fishtree_phylo)
    ),

    ## read trait file, combine with *first* tree block phylo
    tar_target(
        ag_compdata_tb,
        get_ag_data(full_ag_data, phylo = treeblock[[1]])
    )
)  ## end data_input_targets

parameter_constraint_targets <- tar_plan(
        ## set constraints on rate equality
    ## these are derived by staring at the results of corHMM::getStateMat4Dat(ag_compdata$data)
    ## along with state names and figuring out which transitions to constrain

    ## NO constraints (all 24 parameters)
    tar_target(
        full_pars,
        list()
    ),
    ## least-constrained model (12 params: 24 → collapse 4 sets of 4 to 1 each
    tar_target(
        pcsc_pars,
        list(c(7, 10, 20, 23),  ## all pc_gain rates
             c(4, 11, 17, 24),  ## all sc_gain rates
             c(2, 5, 15, 18),   ## all pc_loss rates
             c(1, 8, 14, 21))   ## all sc_loss rates
    ),
    ## add constraints: ag gain/loss depends only on pc
    tar_target(
        pc_pars,
        c(pcsc_pars,
          list(c(13, 16), ## ag_gain with pc==0 (sc==0 or 1)
               c(19, 22), ## ag_gain with pc==1 (sc==0 or 1)
               c(3,   6), ## ag_loss with pc==0
               c(9,  12)) ## ag_gain with pc==1
          )
    ),
    ## add constraints: ag gain/loss depends only on sc
    tar_target(
        sc_pars,
        c(pcsc_pars,
          list(c(13, 19), ## ag_gain with sc==0
               c(16, 22), ## ag_gain with sc==1
               c(3,   9), ## ag_loss with sc==0
               c(6,  12)) ## ag_gain with sc==1
          )
    ),

    ## add constraints: ag gain/loss is independent of sc, pc
    tar_target(
        indep_pars,
        c(pcsc_pars,
          list(c(13, 16, 19, 22), ## ag_gain
               c(3, 6,  9, 12) ## ag_loss
               )
          )
    )
)

##
list(data_input_targets,
     parameter_constraint_targets,

    ## construct state matrices/indices for all models
    tar_map(
        values = tibble(
            ## general trick: pass 'nm' as a column of the values
            ## so we get ag_statemat_{nm}
            nm = c("full", "pcsc", "pc", "sc", "indep"),
            eqstatepars =
              rlang::syms(c("full_pars",
                            "pcsc_pars",
                            "pc_pars",
                            "sc_pars",
                            "indep_pars"))),
        names = nm,
        tar_target(ag_statemat, {
          ag_smdat <- corHMM::getStateMat4Dat(ag_compdata$data)
          equateStateMatPars(ag_smdat$rate.mat, eqstatepars)
        })
    ),

    ## define root state
    tar_target(
        root.p, {
          r <- rep(0, 8)
          r[2] <- 1 ## ag0_care0_spawning1: no ag, no parental care, yes group spawning
          r
        }),

    ## define bounds for corHMM fits
    tar_target(
        ag_corhmm_bounds,
        c(lower = 0.1,      ## 0.1 transitions per tree
          upper = 100 * ape::Ntip(ag_compdata$phy)) ## 100 transitions per species
    ),

    ## fit corHMM models for all sets of constraints
    tar_map(
        values = tibble(
            nm = c("full", "pcsc","pc","sc","indep"),
            statemat = rlang::syms(c("ag_statemat_full",
                                     "ag_statemat_pcsc",
                                     "ag_statemat_pc",
                                     "ag_statemat_sc",
                                     "ag_statemat_indep"))),
        names = nm,
        tar_target(
            ag_model,
            augment_model(
                corHMM(phy = ag_compdata$phy,
                       data = ag_compdata$data,
                       rate.cat = 1,
                       rate.mat = statemat,
                       root.p = root.p,
                       lower = ag_corhmm_bounds[["lower"]],
                       upper = ag_corhmm_bounds[["upper"]]
                       )
                ))
    ),

    tar_target(
        ## FIXME: DRY (via tar_map) and/or don't bother with 'fishphylo' fit?
        ag_model_tb, {
            augment_model(
                ## FIXME: quietly?
                ## (drop node labels for one part)
                corHMM(phy = ag_compdata_tb$phy,
                       data = ag_compdata_tb$data,
                       rate.cat = 1,
                       rate.mat = ag_statemat_pcsc,
                       root.p = root.p,
                       lower = 0.1,                             ## 0.1 transitions per tree
                       upper = 100 * ape::Ntip(ag_compdata$phy) ## 100 transitions per species
                       )
                )
        }),
    tar_target(ag_model_pcsc_prior,
    {
      nll <- make_nllfun(ag_model_pcsc)
      p <- coef(ag_model_pcsc)
      cc <- corhmm_logpostfun
      p[] <- 0
      parnames(cc) <- names(p)
      mle2(cc,
           start = p,
           trace = TRUE,
           vecpar = TRUE,
           data = list(nllfun = nll,
                       negative = TRUE,
                       ## sum(edge length) scaled to 1
                       lb = log(1),
                       ub = log(10 * ape::Ntip(ag_compdata$phy)),
                       gainloss_pairs = gainloss_priors$pairs,
                       lb_gainloss = gainloss_priors$lb,
                       ub_gainloss = gainloss_priors$ub),
           control = list(maxit = 1e4, trace = 1),
           method = "BFGS"
           )
    }),
    tar_target(comp_ci,
    { list(
          wald = tidy(ag_model_pcsc, conf.int = TRUE),
          ## profile = tidy(ag_model_pcsc, conf.int = TRUE,
          ## conf.method = "profile", profile = ag_profile0),
          mcmc = tidy(ag_mcmc0, robust = TRUE, conf.int = TRUE)) %>%
        bind_rows(.id = "method")  %>%
        rename(lwr = "conf.low", upr = "conf.high")
    }
    ),
    tar_target(
        gainloss_priors,
        list(pairs = list(c(4,1), c(6,2), c(9,3), c(10,5), c(11, 7), c(12, 8)),
             ##        sc     pc      ag ....
             ub = log(c(10,    5 ,  rep(10, 4))),
             lb = log(c(0.2,    0.1,  rep(1e-3, 4))))
    ),
    tar_target(
        contrast_mat, {
          cmat <- (read_csv("contr.csv", col_types = cols())
            ## arrange in same order as ag_mcmc1 columns!
            %>% mutate(across(parname, ~factor(., levels=colnames(ag_mcmc1))))
            %>% arrange(parname)
            %>% dplyr::select(-parname)  ## drop row name so we have a pure-numeric matrix
            %>% as.matrix()
          )
          rownames(cmat) <- colnames(ag_mcmc1)
          cmat
        }
    ),
    tar_map(
        values = tibble(mcmc = rlang::syms(c("ag_mcmc0",
                                             "ag_mcmc_tb",
                                             "ag_priorsamp"))),
        tar_target(contr_long,
           ((as.mcmc(mcmc) %*% contrast_mat)
               %>% as_tibble()
               %>% pivot_longer(everything(), names_to = "contrast")
               %>% separate(contrast, into=c("contrast", "rate"))
           ),
        )
    ),
    tar_target(
        states_df, {
          sm <- with(ag_model_pcsc,
                     makeSimmap(phy, data, solution, rate.cat, nSim = 100, nCores = 5))
          purrr::map_dfr(sm, ~ get_state_occ_prop(.[["maps"]])) %>% setNames(state_names(ag_compdata$data[,-1]))
        }),
    tar_target(
        all_ci, {
          t_list <- list(
              wald = tidy(ag_model_pcsc, conf.int = TRUE),
              wald_prior = tidy(ag_model_pcsc_prior, conf.int = TRUE,
                                conf.method = "quad"),
              ## profile = tidy(ag_model_pcsc, conf.int = TRUE,
              ## conf.method = "profile", profile = ag_profile0),
              mcmc = tidy(ag_mcmc0, robust = TRUE, conf.int = TRUE),
              mcmc_treeblock = tidy(ag_mcmc_tb, robust = TRUE, conf.int = TRUE)
          )
          (bind_rows(t_list, .id = "method")
            %>% rename(lwr = "conf.low", upr = "conf.high")
          )
        }),
    tar_target(
        prior_ci,
            (as.mcmc(ag_priorsamp)
                %>% apply(MARGIN = 2, quantile, c(0.025, 0.975))
                %>% t()
                %>% as.data.frame()
                %>% setNames(c("lwr", "upr"))
                %>% tibble::rownames_to_column("term")
            )
    ),
    tar_target(ag_mcmc0,
               corhmm_mcmc(ag_model_pcsc,
                           p_args=list(nllfun = make_nllfun(ag_model_pcsc),
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata$phy)),
                                       gainloss_pairs = gainloss_priors$pairs,
                                       lb_gainloss = gainloss_priors$lb,
                                       ub_gainloss = gainloss_priors$ub),
                           n_cores = 8,
                           n_chains = 8,
                           n_burnin = 4000,
                           n_iter = 84000,
                           n_thin = 20,
                           seed = 101),
               cue = run_mcmc()
               ),
    ## DRY: map with previous rule
    tar_target(ag_mcmc_full,
               corhmm_mcmc(ag_model_full,
                           p_args=list(nllfun = make_nllfun(ag_model_pcsc),
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata$phy)),
                                       gainloss_pairs = gainloss_priors$pairs,
                                       lb_gainloss = gainloss_priors$lb,
                                       ub_gainloss = gainloss_priors$ub),
                           n_cores = 8,
                           n_chains = 8,
                           n_burnin = 4000,
                           n_iter = 84000,
                           n_thin = 20,
                           seed = 101),
               cue = run_mcmc()
               ),
    tar_target(ag_mcmc_tb,
               corhmm_mcmc(ag_model_tb,
                           p_args=list(nllfun = make_nllfun(ag_model_tb, treeblock = treeblock),
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata_tb$phy)),
                                       gainloss_pairs = gainloss_priors$pairs,
                                       lb_gainloss = gainloss_priors$lb,
                                       ub_gainloss = gainloss_priors$ub),
                           n_cores = 8,
                           n_chains = 8,
                           n_burnin =  4000,
                           n_iter =  84000,
                           n_thin = 20,
                           seed = 101),
               cue = run_mcmc()
               ),
    tar_target(ag_priorsamp,
               corhmm_mcmc(ag_model_pcsc,
                           p_args=list(nllfun = function(x) 1,
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata$phy)),
                                       gainloss_pairs = gainloss_priors$pairs,
                                       lb_gainloss = gainloss_priors$lb,
                                       ub_gainloss = gainloss_priors$ub),
                           n_burnin = 4000,
                           n_iter = 84000,
                           n_thin = 20,
                           seed = 101)
               ),
    ## what is this actually doing? doesn't save the objects
    tar_map(
        values = tibble(mcmc = rlang::syms(c("ag_mcmc0", "ag_mcmc_tb"))),
        tar_target(traceplot, lattice::xyplot(mcmc, aspect="fill", layout=c(2,6)))
    ),
    tar_map(values = tibble(mcmc = rlang::syms("ag_mcmc0", "ag_mcmc_tb", "ag_mcmc_full"),
                            nm = c("0", "tb", "full")),
            names = nm,
            tar_target(pairsplots,
                       mk_mcmcpairsplot(mcmc, sprintf("mcmc_pairs_%s.png", nm)))
            ),
    tar_map(
        values = tibble(mcmc = rlang::syms(c("ag_mcmc0", "ag_mcmc_tb")),
                       fn = "pairs_ag_"),
        tar_target(pairsplot, lattice::xyplot(mcmc, aspect="fill", layout=c(2,6)))
    ),
    tar_target(ag_mcmc1,
               as.mcmc(ag_mcmc0)
               ),
    ## SKIP for now
    ## tar_target(ag_profile0,
    ##            profile(ag_model_pcsc,
    ##                    n_cores = 12,
    ##                    trace = TRUE,
    ##                    alpha=0.05) ## less extreme than default (alpha=0.01)
    ##            ),
    tar_render(ag_rmd,
               "ag_model.rmd"
               ),
    tar_render(ag_bayes_rmd,
               "ag_bayesdiag.rmd"
               )
    ## clean up/rescue?
    ## tar_render(ag_model_tech_html,
    ## "ag_model_tech.rmd"
    ## )
)
