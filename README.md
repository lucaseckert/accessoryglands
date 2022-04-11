# accessoryglands

Phylogenetic comparative analysis on the evolution of reproductive accessory glands in fishes.

## workflow

- implemented in the [targets package](https://books.ropensci.org/targets/)
- for different cases (data, settings, etc.), we need to be able to:
   - fit a basic `corHMM` model
   - possibly: compute Wald and/or likelihood profile confidence intervals
   - specify priors and do MCMC computation
   - plot/interpret results

## getting started

```r
source("R/utils.R")
install_pkgs()
```
Files:

- `accessTree.R`
- `binaryTraitData.csv`
- `cleanTraitData.csv`
- `treeBlock.rds`
- `treeSingle.rds`


## making stuff

To run single-threaded:

```r
library(targets)
tar_make()
```

Or multi-threaded:

```r
library(targets)
library(future)
plan("multicore")
tar_make_future()
```

To load existing products (e.g. for debugging):

```r
tar_load(everything())
```
google drive [link](https://drive.google.com/drive/folders/1S5KwLDQavshwS8i0e9_g1jRiVUw8rnLO?usp=sharing)
