#### Artiodactyl BT ####

## packages
library(btw)

## data
dat <- read.table("bayestraits/Artiodactyl.txt")
tt <- ape::read.nexus("bayestraits/Artiodactyl.trees") ## 500 trees

## trying it the normal way
command_vec<-c("1","2","PriorAll exp 10")
r<-bayestraits(data = dat, tree = tt, commands = command_vec) ## cant find BT

## trying it with our functions
## command function
bt_command <- function(prior = NULL,
                       iterations = 51e4,
                       burnin = 1e4,
                       sample = 1,
                       seed = 101,
                       cores = getOption("bayestraits.cores", 1)) {
  cvec <- c("1", ## MultiState
            "2", ## MCMC
  )
  if (!is.null(prior)) cvec <- c(cvec, prior)
  cvec <-  c(cvec,
             sprintf("Iterations %d", iterations),
             sprintf("Burnin %d", burnin),
             sprintf("Seed %d", seed),
             sprintf("Sample %d", sample)  ## thinning
  )
  if (!is.null(cores)) cvec <- c(cvec, sprintf("Cores %d", cores))
  return(cvec)
}

## run function
bt_run <- function(data = NULL, trees = NULL, prior = NULL, dir = "bayestraits",
                   fn = "", chains = 4, seed0 = 100, verbose = FALSE,
                   ...) {
  all_res <- lapply(seq(chains),
                    function(i) {
                      if (verbose) cat("chain ", i, "\n")
                      command <- bt_command(prior = prior, ...,
                                            seed = seed0 + i)
                      tt <- system.time(
                        res <- bayestraits(data, trees, command)
                      )
                      attr(res, "time") <- tt
                      return(res)
                    })
  results <- list(Log = list(options = all_res[[1]]$Log$options, ## same opts for all chains
                             results = do.call(rbind, lapply(seq_along(all_res), function(i) data.frame(chain = i, all_res[[i]]$Log$results)))
  ))
  results[["time"]] <- lapply(all_res, function(x) attr(x, "time"))
  for (c in c("Schedule", "Stones", "AncStates", "OutputTrees")) {
    results[[c]] <- lapply(all_res, function(x) x[[c]])
  }
  if (nchar(fn) > 0) saveRDS(results, file = file.path(dir, fn))
  return(results)
}

## run
r <- bt_run(prior = "PriorAll exp 10", data = dat, trees = tt, verbose = TRUE) ## error
