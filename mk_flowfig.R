source("R/pkg_install.R")
load_pkgs()
library(targets)
library(colorspace)
tar_load(ag_model_pcsc)

library(diagram)
source("R/plotmat.R") ## hacked version for nudging arrows
source("R/utils.R")

## devs <- c("pdf", "svg", "png", "tikz")
devs <- "tikz"  ## all we really need
for (o in devs) {
    f <- get(o)
    f(sprintf("flowfig.%s", o)); mk_flowfig(tikz = (o=="tikz")); dev.off()
    f(sprintf("flowfig2.%s", o)); mk_flowfig(with_labs = TRUE,
                                             tikz = (o == "tikz")); dev.off()
}
for (f in c("flowfig", "flowfig2")) {
    system(sprintf("pdflatex %s.tikz; convert -density 300 %s.pdf %s.png",
                   f, f, f))
}

## https://tex.stackexchange.com/questions/123924/indexed-letters-inside-circles
