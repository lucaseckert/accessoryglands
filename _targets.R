library("targets")
library("tarchetypes")
library("tidyverse")
source("R/utils.R")
source("R/pkg_install.R")
source("R/functions.R")
source("R/mcmc.R")
options(tidyverse.quiet = TRUE)
grDevices::X11.options(type = "cairo")
options(bitmapType = "cairo")
tar_option_set(packages = all_pkgs)

## prevent parallelization from getting out of hand
Sys.setenv(OPENBLAS_NUM_THREADS="1")
Sys.setenv(OMP_NUM_THREADS="1")

## this flag enables/disables the really slow steps (MCMC sampling, MCMC pairs plots)
redo_slow <- FALSE
run_slow <- function() {
    tar_cue(if (redo_slow) "thorough" else "never")
}

corhmm_models <- c("full", "pcsc", "pc", "sc", "indep")
mcmc_runs <- c("0", "tb", "full", "tb_nogainloss")
mcmc_runs_12 <- c("0", "tb", "tb_nogainloss")  ## 12-parameter models only

## TO DEBUG: set tar_option_set(debug = "target_name"); tar_make(callr_function = NULL), e.g.
## tar_option_set(debug = "ag_pcsc_pars")
## tar_option_set(debug = "ag_compdata")

## to play around interactively:
## library(targets); for (f in c("utils", "functions", "mcmc")) source(sprintf("R/%s.R", f))
## tar_load(everything())

## tar_plan is from tarchetypes::tar_plan()

## rules for importing data
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

## rules for setting up parameter constraints
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

## main loop
list(data_input_targets,
     parameter_constraint_targets,

    ## construct state matrices/indices for all models
    tar_map(
        values = tibble(
            ## general trick: pass 'nm' as a column of the values
            ## so we get ag_statemat_{nm}
            nm = corhmm_models,
            eqstatepars =
                rlang::syms(glue::glue("{corhmm_models}_pars"))
        ),
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

    ## fit corHMM models for all sets of constraints ('fishphylo' phylogeny only)
    tar_map(
        values = tibble(
            nm = corhmm_models,
            statemat = rlang::syms(glue::glue("ag_statemat_{corhmm_models}"))
            ),
        names = nm,
        tar_target(
            ag_model,
            augment_model(
                corHMM(phy = ag_compdata$phy,
                       data = ag_compdata$data,
                       rate.cat = 1,
                       rate.mat = statemat,
                       root.p = root.p,
                       lower.bound = ag_corhmm_bounds[["lower"]],
                       upper.bound = ag_corhmm_bounds[["upper"]]
                       )
                ))
    ),

    ## fit corHMM models for all sets of constraints (using tree block)
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

    ## fit corHMM model with priors (MAP estimation)
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

    ## fit additive model by corHMM
    tar_target(ag_model_pcsc_add,
               fit_contrast.corhmm(ag_model_pcsc,
                                   contrast_mat_inv,
                                   ## zero out interactions
                                   fixed_vals = c(pcxsc_loss = 0, pcxsc_gain = 0),
                                   optControl = list(maxit = 20000))
               ),

    ## compute confidence intervals
    tar_target(comp_ci,
    { list(
          wald = tidy(ag_model_pcsc, conf.int = TRUE),
          ## profile = tidy(ag_model_pcsc, conf.int = TRUE,
          ## conf.method = "profile", profile = ag_profile0),
          mcmc = tidy(ag_mcmc_0, robust = TRUE, conf.int = TRUE)) %>%
        bind_rows(.id = "method")  %>%
        rename(lwr = "conf.low", upr = "conf.high")
    }
    ),

    ## define upper bound for gain/loss ratio
    tar_target(
        maxgainloss,
        log(c(sc = 10, pc = 5, 10))
    ),

    ## define lower bound for gain/loss ratio
    tar_target(
        mingainloss,
        log(c(sc = 1/5, pc = 1/10, 1/1000))
    ),

    ## set up gain/loss priors
    tar_target(
        gainloss_priors,
        list(pairs = gainloss_ind_prs(ag_model_pcsc),
             ##        sc     pc   ag ....
             ub = rep(maxgainloss, c(1, 1, 4)),
             lb = rep(mingainloss, c(1, 1, 4)))
    ),
    ## pattern of focal gain/loss trait by elements of full
    ##  contrast matrix (can be automated??)
    tar_target(
        full_pattern,
        c("sc", "pc", "ag", "pc",
          "ag", "sc", "ag", "ag",
          "sc", "pc", "pc", "sc")
    ),

    ## gain/loss parameters for full (24-parameter) model
    tar_target(
        gainloss_priors_full,
        list(pairs = gainloss_ind_prs(ag_model_full),
             ub = maxgainloss[full_pattern],
             lb = mingainloss[full_pattern])
    ),

    ## define column (contrast) order for 12- and 24-parameter models
    tar_map(
        values = tibble(mcmc = rlang::syms(c("ag_mcmc_0", "ag_mcmc_full")),
                        nm = c("0", "full")),
        names = nm,
        tar_target(col_order,
                   colnames(mcmc[[1]]))
    ),

    ## read in contrast matrices and enforce correct ordering
    tar_map(
        values = tibble(fn = c("contr.csv", "contr_full.csv", "contr_invertible.csv"),
                        col_order = rlang::syms(c("col_order_0", "col_order_full", "col_order_0")),
                        nm = c("0", "full", "inv")
                       ),
        names = nm,
        tar_target(contrast_mat,
        {
            cmat0 <- read_csv(fn, col_types = cols())
            cmat <- (cmat0
                ## arrange in same order as ag_mcmc1 columns!
                %>% mutate(across(parname, ~factor(., levels = col_order)))
                %>% arrange(parname)
                %>% dplyr::select(-parname)  ## drop row name so we have a pure-numeric matrix
                %>% as.matrix()
            )
            rownames(cmat) <- col_order
            cmat
        }
        )
    ),

    ## compute contrasts for 12-parameter fits (fishphylo, treeblock,
    ##  treeblock without gain/loss priors)
    tar_map(
        values = tibble(mcmc = rlang::syms(glue::glue("ag_mcmc_{mcmc_runs_12}"))),
        tar_target(contr_long,
           ((as.mcmc(mcmc) %*% contrast_mat_0)
               %>% as_tibble()
               %>% pivot_longer(everything(), names_to = "contrast")
               %>% separate(contrast, into=c("contrast", "rate"))
           ),
        )
    ),

    ## stochastic character mapping
    tar_target(
        states_df, {
          sm <- with(ag_model_pcsc,
                     makeSimmap(phy, data, solution, rate.cat, nSim = 100, nCores = 5))
          purrr::map_dfr(sm, ~ get_state_occ_prop(.[["maps"]])) %>% setNames(state_names(ag_compdata$data[,-1]))
        }),
    tar_target(
        mod_list,
        (tibble::lst(ag_model_pcsc, ag_model_pcsc_prior, ag_mcmc_0, ag_mcmc_tb, ag_mcmc_tb_nogainloss,
                     ag_priorsamp)
            %>% setNames(gsub("ag_", "", names(.)))
        )
    ),

    ## collect confidence intervals
    tar_target(
        all_ci,
        purrr::map_dfr(mod_list, my_tidy, .id = "method")
    ),

    ## get contrast CIs for 24-parameter model
    tar_target(
        full_contr_ci,
        my_tidy(ag_mcmc_full, contrast_mat = contrast_mat_full) %>% mutate(method = "full", .before = 1)
    ),

    ## collect contrast CIs from all models
    tar_target(
        all_contr_ci,
        (purrr::map_dfr(mod_list, my_tidy, .id = "method",
                        contrast_mat = contrast_mat_inv)
            %>% bind_rows(full_contr_ci)
            %>% bind_rows((tidy(ag_model_pcsc_add, conf.int = TRUE)
                %>% mutate(method = "model_pcsc_add")
                %>% rename(lwr = "conf.low", upr = "conf.high")))
            ## FIXME: gsub("model", "corhmm" OR "mle", method) ...
        )
    ),

    ## run 12-parameter model (SLOW)
    tar_target(ag_mcmc_0,
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
                           n_thin = 10,
                           seed = 101),
               cue = run_slow()
               ),

    ## run 24-parameter model
    ## DRY: map with previous rule
    tar_target(ag_mcmc_full,
    {
        ## https://groups.google.com/g/openblas-users/c/W6ehBvPsKTw
        corhmm_mcmc(ag_model_full,
                           p_args=list(nllfun = make_nllfun(ag_model_pcsc),
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata$phy)),
                                       gainloss_pairs = gainloss_priors_full$pairs,
                                       lb_gainloss = gainloss_priors_full$lb,
                                       ub_gainloss = gainloss_priors_full$ub),
                           n_cores = 8,
                           n_chains = 8,
                           n_burnin = 8000,
                           n_iter = 144000,
                           n_thin = 10,
                           seed = 101)
               },
                   cue = run_slow()

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
                           n_thin = 10,
                           seed = 101),
               cue = run_slow()
               ),
    tar_target(ag_mcmc_tb_nogainloss,
               corhmm_mcmc(ag_model_tb,
                           p_args=list(nllfun = make_nllfun(ag_model_tb, treeblock = treeblock),
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata_tb$phy))
                                       ),
                           n_cores = 8,
                           n_chains = 8,
                           n_burnin =  4000,
                           n_iter =  84000,
                           n_thin = 10,
                           seed = 101),
               cue = run_slow()
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
                           n_thin = 10,
                           seed = 101)
               ),
    tar_map(
        values = tibble(mcmc = rlang::syms(glue::glue("ag_mcmc_{mcmc_runs}")),
                        nm = mcmc_runs),
        names = nm,
        tar_target(traceplot, lattice::xyplot(mcmc, aspect="fill", layout=c(2,6)))
    ),
    ## split traceplot in two pieces
    tar_target(traceplot_full1,
               lattice::xyplot(ag_mcmc_full[,1:12], aspect = "fill", layout = c(2,6))
    ),
    tar_target(traceplot_full2,
           lattice::xyplot(ag_mcmc_full[,13:24], aspect = "fill", layout = c(2,6))
    ),
    tar_map(values = tibble(mcmc = rlang::syms(glue::glue("ag_mcmc_{mcmc_runs}")),
                            nm = mcmc_runs),
            names = nm,
            tar_target(mc_pairsplots,
                       mk_mcmcpairsplot(mcmc, fn = sprintf("pix/mcmc_pairs_%s.png", nm)),
                       format = "file",
                       cue = run_slow()
                       )
            ),
    ## tar_map(
    ##     values = tibble(mcmc = rlang::syms(c("ag_mcmc_0", "ag_mcmc_tb")),
    ##                    fn = "pairs_ag_"),
    ##     tar_target(pairsplot, lattice::xyplot(mcmc, aspect="fill", layout=c(2,6)))
    ## ),
    tar_target(ag_mcmc1,
               as.mcmc(ag_mcmc_0)
               ),
    ## SKIP for now
    ## tar_target(ag_profile0,
    ##            profile(ag_model_pcsc,
    ##                    n_cores = 12,
    ##                    trace = TRUE,
    ##                    alpha=0.05) ## less extreme than default (alpha=0.01)
    ##            ),
    ## old(ish) technical model info
    ## tar_render(ag_old_rmd, "ag_model.rmd"),
    ## Bayesian diagnostics (roll into/include in supplementary material?)
    ## tar_render(ag_bayesdiag_html, "ag_bayesdiag.rmd"),
    ## technical note (audience: technical users/computational folks)
    ## tar_render(ag_tech_html, "ag_tech.rmd") ## ,
    ## supplementary material (audience: general, stats enthusiasts)
    tar_render(ag_supp_html, "ag_supp.rmd")
    ## tar_render(ag_supp_docx, "ag_supp.rmd", output_format = "word_document")
)
