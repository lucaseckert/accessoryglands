library("targets")
source("R/utils.R")
source("R/functions.R")
source("R/mcmc.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = pkgList)
## TO DEBUG: set tar_option_set(debug = "target_name"); tar_make(callr_function = NULL)
## tar_option_set(debug = "ag_compdata")
list(
    tar_target(
        ## trait data
        ag_binary_trait_file,
        "data/binaryTraitData.csv",
        format = "file"
    ),
    tar_target(
        treeblock_file,
        "data/treeBlock.rds",
        format = "file"
    ),
    tar_target(
        full_ag_data,
        read_csv(ag_binary_trait_file, col_types = cols())
    ),
    tar_target(
        treeblock,
        readRDS("data/treeBlock.rds") %>% purrr::map(scale_phylo)
    ),
    tar_target(
        fishtree_phylo,
        suppressWarnings(
            fishtree_phylogeny(full_ag_data$species)) %>% scale_phylo()
    ),
    tar_target(
        ## read trait file, grab phylo data from fishtree, combine/trim
        ag_compdata,
        get_ag_data(full_ag_data, phylo = fishtree_phylo)
    ),
    tar_target(
        ag_compdata_treeblock,
        get_ag_data(full_ag_data, phylo = treeblock[[1]])
    ),
    tar_target(
        ## state matrix for AG problem (12 compartments)
        ag_statemat1, {
          ## debugging
          ## data <- ag_compdata$data
          ## with(data, table(ag, care, spawning))
          with(ag_compdata, {
              ag_smdat <- corHMM::getStateMat4Dat(data)
              ## FIXME: name/number match here?
              pars2equal <- list(c(7, 10, 20, 23), c(4, 11, 17, 24), c(2, 5, 15, 18), c(1, 8, 14, 21))
              equateStateMatPars(ag_smdat$rate.mat, pars2equal)
          }
          )
        }),
    tar_target(
        root.p,
        {
          r <- rep(0, 8)
          r[2] <- 1 ## ag0_care0_spawning1: no ag, no parental care, yes group spawning
          r
        }),
    tar_target(
        ag_model0, {
            augment_model(
                ## FIXME: quietly?
                ## (drop node labels for one part)
                corHMM(phy = ag_compdata$phy,
                       data = ag_compdata$data,
                       rate.cat = 1,
                       rate.mat = ag_statemat1,
                       root.p = root.p,
                       lower = 0.1,                             ## 0.1 transitions per tree
                       upper = 100 * ape::Ntip(ag_compdata$phy) ## 100 transitions per species
                       )
                )
        }),
    tar_target(
        ## FIXME: DRY
        ag_model_treeblock, {
            augment_model(
                ## FIXME: quietly?
                ## (drop node labels for one part)
                corHMM(phy = ag_compdata_treeblock$phy,
                       data = ag_compdata_treeblock$data,
                       rate.cat = 1,
                       rate.mat = ag_statemat1,
                       root.p = root.p,
                       lower = 0.1,                             ## 0.1 transitions per tree
                       upper = 100 * ape::Ntip(ag_compdata$phy) ## 100 transitions per species
                       )
                )
        }),
    tar_target(
        gainloss_priors,
        list(pairs = list(c(4,1), c(6,2), c(9,3), c(10,5), c(11, 7), c(12, 8)),
             ##        sc     pc      ag ....
             ub = log(c(10,    5 ,  rep(10, 4))),
             lb = log(c(5,    0.1,  rep(1e-3, 4))))
    ),
    tar_target(ag_mcmc0,
               corhmm_mcmc(ag_model0,
                           p_args=list(nllfun = make_nllfun(ag_model0),
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
                           seed = 101)
               ),
    tar_target(ag_mcmc_treeblock,
               corhmm_mcmc(ag_model_treeblock,
                           p_args=list(nllfun = make_nllfun(ag_model_treeblock, treeblock = treeblock),
                                       ## sum(edge length) scaled to 1
                                       lb = log(1),
                                       ub = log(10 * ape::Ntip(ag_compdata_treeblock$phy)),
                                       gainloss_pairs = gainloss_priors$pairs,
                                       lb_gainloss = gainloss_priors$lb,
                                       ub_gainloss = gainloss_priors$ub),
                           n_cores = 8,
                           n_chains = 8,
                           n_burnin =  4000,
                           n_iter =  84000,
                           n_thin = 20,
                           seed = 101)
               ),
    tar_target(ag_mcmc1,
               as.mcmc(ag_mcmc0)
               ),
    ## SKIP for now
    ## tar_target(ag_profile0,
    ##            profile(ag_model0,
    ##                    n_cores = 12,
    ##                    trace = TRUE,
    ##                    alpha=0.05) ## less extreme than default (alpha=0.01)
    ##            ),
    tar_target(ag_rmd,
               format = "file",
               rmarkdown::render("ag_model.rmd")
               )
)

