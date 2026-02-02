## only package needed before starting this is remotes

CRAN_core_pkgs <- c("abind", "bbmle", "bookdown", "broom.mixed", "caper", "coda", 
                    "colorspace", "corHMM", "cowplot", "diagram", "diversitree", 
                    "emdbook", "fishtree", "GGally", "ggmosaic", "ggnewscale", "ggthemes", 
                    "glue", "hues", "igraph", "latticeExtra", "Matrix", 
                    "nloptr", "numDeriv", "patchwork", "phangorn", "phytools", "ramcmc",
                    "remotes", "tarchetypes", "targets", "tidyverse", "visNetwork")

CRAN_pix_pkgs <- c("tikzDevice", "tinytex")

CRAN_pkgs <- c(CRAN_core_pkgs, CRAN_pix_pkgs)

## packages to install from GitHub (username, reponame)
GH_core_pkgs <- list(c("bbolker","corHMM")        ## hacked/BMB version
                   , c("bbolker","btw")           ## BayesTree interface (hacked/BMB version)
                   , c("easystats", "bayestestR") ## need devel version
                )

GH_pix_pkgs <- list(c("YuLab-SMU","ggtree"))

GH_pkgs <- c(GH_core_pkgs, GH_pix_pkgs)
GH_pkg_nms <- sapply(GH_pkgs, function(x) x[2])

all_pkgs <- c(CRAN_pkgs, GH_pkg_nms)

## redundant with `tar_option_set(packages = pkgList)` ?
load_pkgs <- function(skip = NULL) {
    pkgs <- setdiff(all_pkgs, skip)
    invisible(lapply(pkgs, library, character.only = TRUE))
}


base_pkgs <- "Matrix"
##' install uninstalled pkgs from pkgList and GH_pkgs
##' FIXME: use renv() ?
install_pkgs <- function() {
    if (!require("remotes")) stop("please install the 'remotes' package")
    ip <- installed.packages()
    to_install <- setdiff(CRAN_pkgs,
                          ## skip installed packages and packages also on GH
                          c(rownames(ip),
                            base_pkgs,
                            GH_pkgs))
    install.packages(to_install)
    ## since install_github checks hashes, don't bother checking whether already installed
    lapply(GH_pkgs,  function(x)
        remotes::install_github(paste(x[[1]], x[[2]], sep = "/")))
    return(invisible(NULL))
}

