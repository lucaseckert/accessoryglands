get_ag_data <- function(full_data, species_name = "species", trait_names = c("ag", "care", "spawning"), phylo) {
  full_data <- full_data[c(species_name, trait_names)]
  ## prepare binary trait data for corHMM missing convention
  full_data <- (full_data
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
