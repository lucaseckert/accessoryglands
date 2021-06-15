# accessoryglands

Phylogenetic comparative analysis on the evolution of reproductive accessory glands in fishes.

## workflow

- implemented in the [targets package](https://books.ropensci.org/targets/)
- for different cases (data, settings, etc.), we need to be able to:
   - fit a basic `corHMM` model
   - possibly: compute Wald and/or likelihood profile confidence intervals
   - specify priors and do MCMC computation
   - plot/interpret results

Files:

- `accessTree.R`
- `binaryTraitData.csv`
- `cleanTraitData.csv`
- `treeBlock.rds`
- `treeSingle.rds`

