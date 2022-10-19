source("R/utils.R")
source("R/mcmc.R")
source("R/functions.R")
load_pkgs()
library(targets)
library(colorspace)
`tar_load(ag_model_pcsc)

library(diagram)
M <- ag_model_pcsc$solution
R <- ag_model_pcsc$args.list$rate
dimnames(R) <- dimnames(M)
R[R==13] <- NA
storage.mode(R) <- "character"
R[!is.na(R)] <- ""

## hack (reverse order)
vals <- expand.grid(sc=0:1, pc=0:1, ag=0:1)[,3:1]
## put nodes on flattened hypercube
##  (ag=0 inner square, ag = 1 outer square

box_mag <- vals$ag+0.5
shift <- c(0.3, 0.5)
xval <- with(vals, (sc-0.5)*box_mag)
yval <- with(vals, (pc-0.5)*box_mag)
pos <- cbind(xval, yval)
pos2 <- (pos-min(pos) + shift[1])/(diff(range(pos))*(1+shift[2]))
col_vec <- rep(NA, 12)
up_down_ind <- c(1,2,4,6)
col_vec[up_down_ind] <- gray(c(0, 0.3, 0.6, 0.9))
col_vec[-up_down_ind] <- qualitative_hcl(8)
C <- c(R)
C <- col_vec[C]
dim(C) <- dim(R)

plotmat(R, pos = pos2, xlim = c(-3,3),
        ## arr.lwd = sqrt(M/50),
        box.type = "ellipse", box.prop = 0.5,
        arr.lcol = C, arr.col = C)

dev.copy2pdf(file="flowfig.pdf")

svg("flowfig.svg")
plotmat(R, pos = pos2, xlim = c(-3,3),
        ## arr.lwd = sqrt(M/50),
        box.type = "ellipse", box.prop = 0.5,
        arr.lcol = C, arr.col = C)
dev.off()
## post-editing ...
