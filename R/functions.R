get_ag_data <- function(full_data, species_name = "species", trait_names = c("ag", "care", "spawning"), phylo) {
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
  cc <- caper::comparative.data(phy = phylo,
                                data = as.data.frame(full_data), ## f'n works REALLY BADLY with tibbles!
                                names.col = "species")
  ## construct new version of data including species name (corHMM format)
  cc$data <- data.frame(species = rownames(cc$data), cc$data)
  rownames(cc$data) <- NULL
  cc$data <- shorten_ag_names(cc$data)
  return(cc)
}

## MINOR change to get_ag_data screws up lots of downstream stuff??
## create a new function to stay safe?
get_ag_data2 <- function(full_data, species_name = "species", trait_names = c("ag", "care", "spawning"), phylo) {
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
  cc <- caper::comparative.data(phy = phylo,
                                data = as.data.frame(full_data), ## f'n works REALLY BADLY with tibbles!
                                names.col = "species")
  ## construct new version of data including species name (corHMM format)
  cc$data <- data.frame(species = rownames(cc$data), cc$data)
  rownames(cc$data) <- NULL
  cc$data <- shorten_ag_names(cc$data)
  return(cc)
}
