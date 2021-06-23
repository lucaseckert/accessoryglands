get_ag_data <- function(file, species_name = "species", trait_names = c("ag", "care", "spawning")) {
  full_data <- (file
    %>% read_csv(col_types = cols())
  )
  full_data <- full_data[c(species_name, trait_names)]
  ## %>% dplyr::select(species, ag, care, spawning)
  ## prepare binary trait data for corHMM missing convention
  full_data <- (full_data
    %>% tidyr::replace_na(list(ag = "?",
                               care = "?",
                               spawning = "?"))
  )
  phylo <-  suppressWarnings(
      fishtree_phylogeny(full_data$species))
  ## request 607, only find 478 spp
  cc <- caper::comparative.data(phy = phylo,
                          data = as.data.frame(full_data), ## f'n works REALLY BADLY with tibbles!
                          names.col = "species")
  ## construct new version of data including species name (corHMM format)
  cc$data0 <- cc$data
  cc$data <- data.frame(species = rownames(cc$data0), cc$data0)
  rownames(cc) <- NULL
  return(cc)
}
