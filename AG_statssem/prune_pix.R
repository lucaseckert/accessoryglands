#' Plot a vertical cladogram and overlay horizontal state-probability boxes
#'
#' Tip rows are drawn just ABOVE the tips.
#' Internal-node rows are drawn just BELOW the node (on the incoming branch).
#'
#' @param tree ape::phylo
#' @param P numeric matrix, nrow = Ntip(tree) + tree$Nnode, ncol = K states.
#'          Row i corresponds to node i in standard phylo numbering:
#'          1:Ntip(tree) = tips; Ntip+1 : Ntip+Nnode = internal nodes.
#' @param digits printed digits for probabilities
#' @param tip_offset_y fraction of y-range to shift tip rows UP
#' @param node_offset_y fraction of y-range to shift internal-node rows DOWN
#' @param box_height_y fraction of y-range for box height
#' @param box_width_x fraction of x-range for cell width (each state cell)
#' @param box_gap_x fraction of x-range between adjacent state cells
#' @param border,fill colors; NA cells are left blank (no text)
#' @param cex text size inside cells
#' @param plot_args list of extra args for ape::plot.phylo
#'        (e.g. edge.width, show.tip.label, etc.)
#' @return (invisible) list with node coordinates and box geometry
plot_pruning_boxes_vertical <- function(
  tree, P,
  digits       = 2,
  tip_offset_y  = 0.02,
  node_offset_y = 0.02,
  box_height_y  = 0.045,
  box_width_x   = NULL,   # if NULL, auto-scale by K
  box_gap_x     = 0.006,
  border        = "black",
  fill          = "white",
  cex           = 0.8,
  plot_args     = list(direction = "upwards",
                       type = "cladogram",
                       label.offset = 0.01,
                       show.tip.label = TRUE)
) {
  stopifnot(inherits(tree, "phylo"))
  if (!is.matrix(P)) P <- as.matrix(P)
  N  <- ape::Ntip(tree)
  Ni <- tree$Nnode
  n_total <- N + Ni
  if (nrow(P) != n_total)
    stop("nrow(P) must equal Ntip(tree) + tree$Nnode")

  K <- ncol(P)

  ## 1) Plot the vertical cladogram
  do.call(ape::plot.phylo, append(list(x = tree), plot_args))
  last <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  xx <- last$xx  # horizontal coordinate
  yy <- last$yy  # vertical coordinate (root near min, tips near max)

  ## 2) Convert relative sizing to device units
  usr   <- par("usr")         # c(x1, x2, y1, y2)
  xspan <- usr[2] - usr[1]
  yspan <- usr[4] - usr[3]

  if (is.null(box_width_x)) {
    # pick a reasonable width per cell and adapt to K
    box_width_x <- (0.12 * xspan) / max(K, 3)
  }
  bw  <- box_width_x
  bh  <- box_height_y * yspan
  gap <- box_gap_x * xspan

  # Helper: draw K horizontally arranged cells centered at y_mid,
  # starting at x_left, all at the same y.
  draw_row <- function(x_left, y_mid, probs) {
    for (j in seq_len(K)) {
      xl <- x_left + (j - 1) * (bw + gap)
      rect(xl, y_mid - bh/2, xl + bw, y_mid + bh/2,
           col = fill, border = border)
      if (!is.na(probs[j])) {
        text(xl + bw/2, y_mid,
             formatC(probs[j], format = "f", digits = digits),
             cex = cex)
      }
    }
  }

  ## 3) Tips: draw rows just ABOVE the tip (positive y offset)
  for (i in seq_len(N)) {
    x_center <- xx[i]
    # Left-align the row symmetrically around the node in x
    row_span <- K * bw + (K - 1) * gap
    x_left   <- x_center - row_span / 2
    y_above  <- yy[i] + tip_offset_y * yspan
    draw_row(x_left, y_above, P[i, ])
  }

  ## 4) Internal nodes: draw rows just BELOW the node (negative y offset)
  for (i in (N + 1):n_total) {
    x_center <- xx[i]
    row_span <- K * bw + (K - 1) * gap
    x_left   <- x_center - row_span / 2
    y_below  <- yy[i] - node_offset_y * yspan
    draw_row(x_left, y_below, P[i, ])
  }

  invisible(list(coords = data.frame(node = 1:n_total, x = xx, y = yy),
                 box = list(cell_width_x = bw, height_y = bh, gap_x = gap, K = K)))
}

## https://blog.phytools.org/2023/03/simple-demonstration-of-felsensteins.html
##' @param q rates
##' @param tree
##' @param x tip states
##' @param model matrix of transition-rate indices
##' optional pi: root priors
pruning <- function(q, tree, x, model=NULL, steps = NULL, ...){
  require("phytools", quietly = TRUE) ## for to.matrix
  pw <- reorder(tree,"postorder")
  k <- length(levels(x))
  if (is.null(model)){
    model <- matrix(1,k,k)
    diag(model) <- 0
  }
  if(hasArg(pi)) pi <- list(...)$pi
  else pi <- rep(1/k,k)
  Q <- matrix(0,k,k)
  Q[] <- c(0,q)[model+1]
  diag(Q) <- -rowSums(Q)
  L <- rbind(to.matrix(x[pw$tip.label],levels(x)),
      matrix(0,tree$Nnode,k,dimnames=
                              list(1:tree$Nnode+Ntip(tree))))
  nn <- unique(pw$edge[,1]) 
  for (i in 1:length(nn)) {
    ee <- which(pw$edge[,1]==nn[i])
    PP <- matrix(NA,length(ee),k)
    for(j in 1:length(ee)) {
      P <- Matrix::expm(Q*pw$edge.length[ee[j]])
      PP[j,] <- drop(as.matrix(P%*%L[pw$edge[ee[j],2],]))
    }
    if (!is.null(steps) && i > steps) return(L)
    L[nn[i],] <- apply(PP, 2, prod)
  }
  if (!is.null(steps)) return(L)
  prob <- log(sum(pi*L[nn[i],]))
  prob
}

library(ape)

# A small 6-tip tree with varied branch lengths
tr  <-  read.tree(text = "((A:1,B:1):0.5,(C:1.5,(D:0.5,E:0.5):2.0):1.0,F:2.5):0.5;")


set.seed(101)
invisible(capture.output(
  tr  <-  rtree(6) |> phytools::force.ultrametric() |>
    reorder("postorder")
))
tr$edge.length <- round(tr$edge.length/min(tr$edge.length))
N   <-  Ntip(tr)
Ni  <-  tr$Nnode
K   <-  3
ntot  <-  N + Ni

# Start with all NA, then fill tips with one-hot encoding of observed states
P  <-  matrix(NA_real_, nrow = ntot, ncol = K)

# Tip states: random
tip_states  <-  sample(0:(K-1), N, replace = TRUE)
for (i in 1:N) P[i, ]  <-  as.numeric(0:(K - 1) == tip_states[i])

pruning(q=rep(1,6), tree = tr, x = setNames(factor(tip_states), tr$tip.label))
pruning(q=rep(1,6), tree = tr, x = setNames(factor(tip_states), tr$tip.label), steps = 1)
pruning(q=rep(1,6), tree = tr, x = setNames(factor(tip_states), tr$tip.label), steps = 2)
pruning(q=rep(1,6), tree = tr, x = setNames(factor(tip_states), tr$tip.label), steps = 3)
pruning(q=rep(1,6), tree = tr, x = setNames(factor(tip_states), tr$tip.label), steps = 4)
  
par(mar = c(2,2,2,2), xpd = NA)
for (i in 0:5) {
  P <- pruning(q=rep(1,6), tree = tr, x = setNames(factor(tip_states), tr$tip.label), steps = i)
  P[rowSums(P)==0,] <- NA
  res  <-  plot_pruning_boxes_vertical(
    tr, P,
    digits = 3,
    box_gap_x = 0,
    box_width_x = 0.3
  )
  
}

nodelabels()
## run the pruning algorithm!

# Partially fill some internal nodes (to mimic "has reached/not yet")
# NOTE: internal nodes are numbered N+1 to N+Ni
P[N + 1, ]  <-  c(0.5, NA, NA)
P[N + 2, ]  <-  c(NA, 1.0, NA)
P[N + 3, ]  <-  c(NA, NA, NA)   # still blank
P[N + 4, ]  <-  c(0.5, NA, NA)

# Plot
title(main = "Pruning algorithm: probabilities filled where available", adj = 0)
