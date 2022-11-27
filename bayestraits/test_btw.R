remotes::install_github("bbolker/btw")
library(btw)
library(targets)
library(tidyverse)

## devtools::load_all("~/R/pkgs/btw")
options(bt_path = "./BayesTraitsV4.0.0-Linux",
        bt_bin = "BayesTraitsV4")
bayestraits(data = primate.discrete1,
	tree = primate.tree1,
	command = c("1", "1", "run"))

tar_load(treeblock)
tbmp <- do.call(c, treeblock)
tbmp <- .compressTipLabel(tbmp)
tar_load(ag_compdata_tb)
dd <- (ag_compdata_tb$data
    %>% mutate(across(c("ag", "pc", "sc"),
                      ~ifelse(. == "?", "-", .)))
)
bayestraits(data  = dd,
            tree = tbmp,
            command = c(
                "3", ## dependent
                "2", ## MCMC
                "ScaleTrees",
                "AddTag Root Erpetoichthys_calabaricus Mugil_liza",
                "Fossil Node01 Root 0",
                "Iterations 100",
                "Run"
                ))


