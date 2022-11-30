## only package needed before starting this is remotes

pkgList <- c("abind", "bbmle", "bookdown", "broom.mixed", "caper", "coda", 
             "colorspace", "corHMM", "cowplot", "diagram", "diagram", "diversitree", 
             "emdbook", "fishtree", "GGally", "ggmosaic", "ggnewscale", "ggthemes", 
             "ggtree", "glue", "hues", "igraph", "latticeExtra", "Matrix", 
             "nloptr", "numDeriv", "patchwork", "phangorn", "phytools", "ramcmc", "remotes", 
             "tarchetypes", "targets", "tidyverse", "tikzDevice", "visNetwork")

## packages to install from GitHub (username, reponame)
GH_pkgs <- list(c("bbolker","corHMM"),    ## hacked/BMB version
                c("YuLab-SMU","ggtree"),
                c("bbolker","btw")       ## BayesTree interface (hacked/BMB version)
                )

base_pkgs <- "Matrix"
##' install uninstalled pkgs from pkgList and GH_pkgs
##' FIXME: use renv() ?
install_pkgs <- function() {
    if (!require("remotes")) stop("please install the 'remotes' package")
    ip <- installed.packages()
    to_install <- setdiff(pkgList,
                          ## skip installed packages and GitHub packages
                          c(rownames(ip),
                            sapply(GH_pkgs, function(x) x[2])))
    install.packages(to_install)
    ## since install_github checks hashes, don't bother checking whether already installed
    lapply(GH_pkgs,  function(x)
        remotes::install_github(paste(x[[1]], x[[2]], sep = "/")))
    return(invisible(NULL))
}

