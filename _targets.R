library("targets")
source("R/utils.R")
source("R/functions.R")
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
        ## read trait file, grab phylo data from fishtree, combine/trim
        ## FIXME: alternative target will be needed if we want to do tree
        ##  blocks (trimming procedure will be different)
        ag_compdata,
        { get_ag_data(
            ag_binary_trait_file
          )
        }
    ),
    tar_target(
        ## state matrix for AG problem (12 compartments)
        ag_statemat1, {
          ## debugging
          ## data <- ag_compdata$data
          ## with(data, table(ag, care, spawning))
          with(ag_compdata, {
            ag_smdat <- corHMM::getStateMat4Dat(data)
            pars2equal <- list(c(7, 10, 20, 23), c(4, 11, 17, 24), c(2, 5, 15, 18), c(1, 8, 14, 21))
            StateMatA_constrained <- equateStateMatPars(ag_smdat$rate.mat, pars2equal)
          }
          )
        }),
    tar_target(
        ag_model0, {
            augment_model(
                corHMM(phy = ag_compdata$phy,
                       data = ag_compdata$data,
                       rate.cat = 1,
                       rate.mat = ag_statemat1)
                )
        }),
    tar_target(
      ag_exp_outfile,
      format = "file",
      {
        fn <-  "cache/ag_fit.rda"
        save(ag_compdata, ag_model0, file = fn)
        fn ## must return fn
      }
    )
)
