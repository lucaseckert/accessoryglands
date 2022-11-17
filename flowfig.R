source("R/utils.R")
source("R/mcmc.R")
source("R/functions.R")
load_pkgs()
library(targets)
library(colorspace)
tar_load(ag_model_pcsc)

library(diagram)
source("R/plotmat.R")
M <- ag_model_pcsc$solution
R <- ag_model_pcsc$args.list$rate
dimnames(R) <- dimnames(M)
R[R==13] <- NA
col_vec <- rep(NA, 12)
up_down_ind <- c(1,2,4,6)
col_vec[up_down_ind] <- gray(c(0, 0.3, 0.6, 0.9))
col_vec[-up_down_ind] <- qualitative_hcl(8)
C <- c(R)
C <- col_vec[C]
dim(C) <- dim(R)

## nudges
nmag <- 0.04

mksymm <- function(X) {
    eqzero <- function(z) !is.na(z) && z == 0
    for (i in 1:nrow(X)) {
        for (j in 1:ncol(X)) {
            if (eqzero(X[i,j]) && !eqzero(X[j,i])) {
                X[i,j] <- X[j,i]
            }
            if (eqzero(X[j,i]) && !eqzero(X[i,j])) {
                X[j,i] <- X[i,j]
            }
        }
    }
    X
}

Nx <- M
Nx[!is.na(Nx)] <- 0
Ny <- Nx

## outer up/down arrows: nudge L/R
Nx[matrix(c("ag1_pc1_sc0", "ag1_pc0_sc0",
            "ag1_pc1_sc1", "ag1_pc0_sc1"),
          byrow = TRUE, ncol = 2)] <- c(-nmag, nmag)
Nx <- mksymm(Nx)
Ny[matrix(c("ag1_pc1_sc0", "ag1_pc1_sc1",
            "ag1_pc0_sc0", "ag1_pc0_sc1"),
          byrow = TRUE, ncol = 2)] <- c(nmag, -nmag)
Ny <- mksymm(Ny)

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

mkplot <- function()
    plotmat(R, pos = pos2, xlim = c(-3,3),
            ## arr.lwd = sqrt(M/50),
            box.type = "ellipse", box.prop = 0.5,
            arr.lcol = C, arr.col = C,
            arr.nudge.x = Nx,
            arr.nudge.y = Ny)

pdf("flowfig.pdf"); mkplot(); dev.off()
svg("flowfig.svg"); mkplot(); dev.off()
png("flowfig.png"); mkplot(); dev.off()
