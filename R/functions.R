get_ag_data <- function(file) {
  full_data <- (file
    %>% read_csv(col_types = cols())
    %>% dplyr::select(species, ag, care, spawning)
    ## prepare binary trait data for corHMM missing convention
    %>% tidyr::replace_na(list(ag = "?",
                               care = "?",
                               spawning = "?"))
  )
  phylo <-  suppressWarnings(
      fishtree_phylogeny(full_data$species))
  ## request 607, only find 478 spp
  caper::comparative.data(phy = phylo,
                          data = as.data.frame(full_data), ## f'n works REALLY BADLY with tibbles!
                          names.col = "species")
}
