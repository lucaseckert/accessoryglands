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
col_vec1 <- rep(NA, 12)
## match by name ...
## up_down_ind == all arrows between equal ag*
ag_labs <- substr(colnames(R), 1, 3)
updown <- outer(ag_labs, ag_labs, "==")
up_down_ind <- sort(unique(R[!is.na(R) & updown])) ## c(1,2,4,6) 
col_vec1[up_down_ind] <- gray(c(0, 0.3, 0.6, 0.9))
col_vec1[-up_down_ind] <- qualitative_hcl(8)
C1 <- c(R)
C1 <- col_vec1[C1] ## map colors to indices
dim(C1) <- dim(R)

## constrained model colour vectors
## pc-only model
col_vec2 <- col_vec1
pc_pairs <- list(c(5, 6), c(7, 8), c(1, 2), c(3, 4))
for (i in seq_along(pc_pairs)) {
    col_vec2[pc_pairs[[i]][2]] <- col_vec2[pc_pairs[[i]][1]]
}
C2 <- c(R)
C2 <- col_vec2[C2] ## map colors to indices
dim(C2) <- dim(R)

col_vec3 <- col_vec1
sc_pairs <- list(c(5, 7), c(6, 8), c(1, 3), c(2, 4))
for (i in seq_along(sc_pairs)) {
    col_vec3[sc_pairs[[i]][2]] <- col_vec3[sc_pairs[[i]][1]]
}
C3 <- c(R)
C3 <- col_vec2[C3] ## map colors to indices
dim(C3) <- dim(R)

                                          

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
## outer L/R arrows: nudge up/down
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

R2 <- R
from <- sprintf("ag0_pc%d_sc%d", rep(0:1,2), rep(0:1, each=2))
to <-   sprintf("ag1_pc%d_sc%d", rep(0:1,2), rep(0:1, each=2))
R2[cbind(from, to)] <- LETTERS[1:4]

mkplot <- function(R, C = C1) {
    plotmat(R, pos = pos2, xlim = c(-3,3),
            ## arr.lwd = sqrt(M/50),
            box.type = "ellipse", box.prop = 0.5,
            arr.lcol = C, arr.col = C,
            arr.nudge.x = Nx,
            arr.nudge.y = Ny)
}

## with labels, for supp/contrast explanation

mkplot2 <- function(R = R2, C = C1) {
    pp <- mkplot(R, C)
    w <- nzchar(na.omit(c(R)))
    ## not quite sure why this fussing is required
    nx <- 0.015*c(0.75,0.75,-0.75,-1)
    ny <- 0.015*c(-0.9,0.9,0.9,0)
    with(pp$arr, mapply(plotrix::draw.circle, TextX[w]+nx, TextY[w]+ny,
                        MoreArgs=list(radius=0.025)))
    invisible(NULL)
}

for (o in c("pdf", "svg", "png")) {
    f <- get(o)
    f(sprintf("flowfig.%s", o)); mkplot(R); dev.off()
    f(sprintf("flowfig2.%s", o)); mkplot2(); dev.off()
}

## https://tex.stackexchange.com/questions/123924/indexed-letters-inside-circles
