source("R/pkg_install.R")
## install_pkgs() if necessary
source("R/utils.R")
load_pkgs()
library(targets)
library(colorspace)
tar_load(ag_model_pcsc)

library(diagram)
source("R/plotmat.R") ## hacked version for nudging arrows

f_args <- list(shadow.size=0, arr.type = "triangle", arr.length = 0.25, arr.width = 0.2, box.col = "cadetblue1", box.lcol = NA)

do.call(mk_flowfig, f_args)

devs <- c("pdf", "svg", "png", "tikz")
## devs <- "tikz"  ## all we really need
for (o in devs) {
    f <- get(o)
    f(sprintf("flowfig.%s", o));
    do.call(mk_flowfig, c(list(tikz = (o=="tikz")), f_args))
    dev.off()
    f(sprintf("flowfig2.%s", o))
    mk_flowfig(with_labs = TRUE,
               do.call(mk_flowfig, c(list(tikz = (o=="tikz"),
                                          with_labs = TRUE),
                                     f_args))); dev.off()
}
for (f in c("flowfig", "flowfig2")) {
    system(sprintf("pdflatex %s.tikz; convert -density 300 %s.pdf %s.png",
                   f, f, f))
}

## https://tex.stackexchange.com/questions/123924/indexed-letters-inside-circles

