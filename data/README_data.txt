DATA FOR ANALYSIS OF ACCESSORY GLAND EVOLUTION

DESCRIPTION OF FILES:

cleanTraitData.csv: trait data for all species in full format
-order: the order the species belongs to
-family: the family the species belongs to
-genus the genus the species belongs to
-species: the species name, corresponding to our phylogeny
-common.name: the common name of the species
-ag.presence: accessory gland presence/absence, 0 if absent, 1 if present
-ag.origin: the structure the gland originates from, S for sperm duct, T for testes, H if hermaphroditic, n/a if no gland is present
-ag.type: the type of accessory gland as determined by the structure it associates with, S for sperm duct gland, M for mesorchial gland, T for testicular gland, D for dorsal accessory gland, n/a if no gland present
-env: the environment the species inhabits, F for freshwater, M for marine
-fertilization: mode of fertilization, E for external, I for internal
-mating: mating system, M for monogamy, G for polygyny, A for polyandry, Y for polygyandry, C for promiscuity, H for hermaphroditic species
-care: form of parental care, P for male, F for female, B for biparental, n/a for no care
-spawning: spawning mode, P for pairs, G for groups, S for alternative reproductive tactics (ARTs)


binaryTraitData.csv: binary version of trait data for analysis
-order: the order the species belongs to
-family: the family the species belongs to
-genus the genus the species belongs to
-species: the species name, corresponding to our phylogeny
-ag: presence/absence of accessory glands, 0 for absent, 1 for present
-mating: NOT USED IN ANALYSIS
-care: form of care, 0 for no care or female care, 1 for male care or biparental care
-spawning: spawning mode, 0 for pair spawning, 1 for groups or ARTs
-pairs.groups: NOT USED IN ANALYSIS
-pairs.arts: NOT USED IN ANALYSIS

dataReferences.csv: references for data for all species
-order: the order the species belongs to
-family: the family the species belongs to
-genus the genus the species belongs to
-species: the species name
-FishTree Name: name corresponding to our phylogeny
-Common Name: common name of the species
-references: all references for that species for all traits in cleanTraitData.csv

referenceList.rtf: full reference list for the sources listed in dataReferences.csv

accessTree.R: code for accessing the tree block and single phylogeny from Fish Tree of Life

treeBlock.rds: R object of tree block (100 phylogenies) for all species in our dataset, generated from accessTree.R

treeSingle.rds: R object of fully resolved phylogeny for species in our dataset, generated from accessTree.R
