source("R/utils.R")
source("R/pkg_install.R")
source("R/functions.R")
source("R/mcmc.R")
library(targets)
library(mirai)
library(parallel)
library(tidyverse)
library(magrittr)  ## need for fancy %>% use
OkIt <- unname(palette.colors(n = 9, palette = "Okabe-Ito"))[-1]

options(tidyverse.quiet = TRUE)
grDevices::X11.options(type = "cairo")
options(bitmapType = "cairo")
tar_option_set(packages = all_pkgs)

library(corHMM)
tar_load(ag_model_pcsc)
tar_load(ag_corhmm_bounds)
tar_load(ag_compdata)
tar_load(gainloss_priors)
tar_load(pcsc_pars)
tar_load(ag_statemat_pcsc)

tar_load(ag_corhmm_bounds)
tar_load(root.p)
tar_load(contrast_mat_inv)
tar_load(full_ag_data)
tar_load(fishtree_phylo)

## refit ag_model_pcsc with priors

fit_model <- function(data) {
    cat("fitting base model\n")
    model0 <- augment_model(
        corHMM(phy = data$phy,
               data = data$data,
               rate.cat = 1,
               rate.mat = ag_statemat_pcsc,
               root.p = root.p,
               lower.bound = ag_corhmm_bounds[["lower"]],
               upper.bound = ag_corhmm_bounds[["upper"]]
               )
    )


    ## re-fit with priors
    nll <- make_nllfun(model0)
    p <- coef(model0)
    cc <- corhmm_logpostfun
    p[] <- 0
    parnames(cc) <- names(p)

    cat("fitting model with priors\n")
    model1 <- mle2(cc,
                   start = p,
                   trace = TRUE,
                   vecpar = TRUE,
                   data = list(nllfun = nll,
                               negative = TRUE,
                               ## sum(edge length) scaled to 1
                               lb = log(1),
                               ub = log(10 * ape::Ntip(ag_compdata$phy)),
                               gainloss_pairs = gainloss_priors$pairs,
                               lb_gainloss = gainloss_priors$lb,
                               ub_gainloss = gainloss_priors$ub),
                   control = list(maxit = 1e4, trace = 1),
                   method = "BFGS"
                   )
    coefs <- my_tidy(model1, contrast_mat = contrast_mat_inv)
    ## now get params with confidence intervals
    return(tibble::lst(model0, model1, coefs))
}

## various subset cases ...

## updated version (for use when switching branches)
get_ag_data <- function(full_data, species_name = "species", trait_names = c("ag", "care", "spawning"), phylo,
                        include_missing = TRUE, subsample = 1.0) {
  full_data <- full_data[c(species_name, trait_names)]
  ## prepare binary trait data for corHMM missing convention
  full_data <- (full_data
    ## make sure data are character ("0"/"1"/NA, not 0/1/NA)
    %>% mutate(across(where(is.numeric), as.character))
    ## corHMM wants NA as ?
    %>% tidyr::replace_na(list(ag = "?",
                               care = "?",
                               spawning = "?"))
  )
  if (!include_missing) {
      full_data <- filter(full_data, !(ag=="?" | care == "?" | spawning == "?"))
  }
  if (subsample < 1) {
      n_full <- nrow(full_data)
      nsamp <- round(n_full*subsample)
      full_data <- dplyr::slice(full_data, sample(n_full, size = nsamp, replace = FALSE))
  }
  cc <- caper::comparative.data(phy = phylo,
                                data = as.data.frame(full_data), ## f'n works REALLY BADLY with tibbles!
                                names.col = "species")
  ## construct new version of data including species name (corHMM format)
  cc$data <- data.frame(species = rownames(cc$data), cc$data)
  rownames(cc$data) <- NULL
  cc$data <- shorten_ag_names(cc$data)
  return(cc)
}


## system.time(
##     fit_model(ag_compdata)
## )

ag_compdata_nomiss <- get_ag_data(full_ag_data, phylo = fishtree_phylo, include_missing = FALSE)

set.seed(101)
ag_subsamp1 <- 
    replicate(5, get_ag_data(full_ag_data, phylo = fishtree_phylo,
                             subsamp = 0.43), simplify = FALSE)
names(ag_subsamp1) <- paste0("subsamp1_0.43_", 1:5)
subframe <- expand.grid(rep = 1:5,
                        s = seq(0.1, 0.9, by = 0.2))


ag_subsamp2 <- list()
for (i in 1:nrow(subframe)) {
    ag_subsamp2[[i]] <-
        get_ag_data(full_ag_data, phylo = fishtree_phylo,
                    subsamp = subframe$s[i])
}
subsamp2_nms <- sprintf("subsamp2_%1.2f_%d", subframe$s, subframe$rep)
names(ag_subsamp2) <- subsamp2_nms

all_data <- c(list(
    orig = ag_compdata,
    nomiss = ag_compdata_nomiss),
    ag_subsamp1,
    ag_subsamp2)

clust <- mirai::make_cluster(n = 12)
pkgs <- c("bbmle", "corHMM", "broom.mixed")
clusterExport(clust, varlist = c("augment_model",
                                 "ag_compdata",
                                 "ag_corhmm_bounds",
                                 "root.p", "ag_statemat_pcsc",
                                 "corhmm_logpostfun", "make_nllfun",
                                 "gainloss_priors", "pkgs",
                                 "state_names",
                                 "par_names",
                                 "my_tidy",
                                 "contrast_mat_inv",
                                 "fix_cnms"
))
clusterEvalQ(clust, sapply(pkgs, library, character.only = TRUE))

all_fits <- parLapply(clust , X = all_data, fun = fit_model)
## names(all_fits)[8:32] <- subsamp2_nms

## save/load loop
saveRDS(all_fits, "subsamp_fits.rds")
all_fits <- readRDS("subsamp_fits.rds")

all_coefs <- map_dfr(all_fits, ~ .[["coefs"]], .id = "sample")

get_subsampsize <- function(s) {
    case_when(s == "orig" ~ 1.0,
              s == "nomiss" ~ 0.431,
              TRUE ~ as.numeric(stringr::str_extract(s,"0.[0-9]{2}")))
}

keep_interesting <- . %>%  filter(stringr::str_detect(term, "(pc|sc|pcxsc)_(loss|gain)$"))

library(ggplot2); theme_set(theme_bw())
ac2 <- all_coefs |>
    keep_interesting() |>
    mutate(s = get_subsampsize(sample),
           batch = stringr::str_extract(sample, "^[[:alpha:]]+")) |>
    arrange(s) |>
    mutate(across(c(sample, term), ~ forcats::fct_inorder(factor(.))))

ggplot(ac2, aes(x=estimate, y = sample,
                colour = factor(s), group = sample)) +
    geom_point() +
    geom_linerange(aes(xmin=lwr, xmax = upr), size = 2) + 
    scale_colour_manual(values = OkIt, guide = guide_legend(reverse = TRUE),
                        name = "subsample proportion") +
    facet_wrap(~term, scale = "free_y") +
    geom_vline(lty=2, xintercept = 0) +
    geom_vline(lty=3, xintercept = log(c(0.1, 10)))

ggsave("subsamples.pdf", width = 16, height = 8)

## comparing runs here and elsewhere
tar_load(all_contr_ci)
all_contr_ci |> filter(method == "model_pcsc_prior") |> keep_interesting() |>
    select(-c(method, std.error, statistic))

ac2 |> filter(sample=="orig") |> keep_interesting() |>
    select(-c(sample, std.error, statistic))


