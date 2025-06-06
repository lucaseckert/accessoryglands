# accessoryglands

Code and data for analyzing the evolution of reproductive accessory glands across fishes. Specifically, we test whether the evolution of accessory glands is driven by other reproductive traits, namely parental care and spawning mode. 

## Caution: project version

There are some glitches somewhere in the `targets` workflow, or the interaction of the `targets` workflow, that we are trying to understand. Some changes between Nov 2022 and March 2025 (which sounds like a long time, but the project was largely inactive during this time) break some of the sensitivity plots of the MLE estimates (e.g. `ag_pcsc_prior` results no longer approximately match the primary MCMC results.

Checking out a 'salvage' version from Nov 30 2022 and rebuilding the supp appears to work. Re-making the `ag_model_*` objects breaks things, for some reason.

We are still trying to diagnose this.

* different versions of data (treeblock vs fishphylo_tree)?
* something about `targets`/`tarchetypes` versions?
* contrast matrix messed up somehow?

In any case, if you want to rebuild the supplementary material, do *not* run `targets::tar_make()`.

## running the project

### getting started

A recent version of R (i.e. 4.2.0 or higher) and updating your packages is strongly recommended (`update.packages()`). You will probably need to have compilation tools installed on your system, for installing R packages from source.

All the code should run, but compiling everything (including PDF version of the supplementary material) needs a few additional tools: a working LaTeX installation, the [optipng](https://optipng.sourceforge.net/) program, ...

```r
if (!require("remotes")) install.packages("remotes")
detach("package:remotes")
source("R/pkg_install.R")
install_pkgs()
```

### running the workflow

These instructions will run all the code and build the supplementary material.

#### within R

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

#### via `make`

- `make all` if you have `make` installed.

To load existing products (e.g. for debugging):

```r
targets::tar_load(everything())
```

To visualize the workflow:

```r
targets::tar_visnetwork()
```

## figures

- figure 1 (flow diagram): `mk_flowfig.R` (needs LaTeX, imagemagick for PDF to PNG conversion)
- figure 2 (phylogeny): used `diversitree` package, `phytools` for stochastic mapping; overlaid in PowerPoint
- figure 3: run `figures.R`

## workflow

- implemented in the [targets package](https://books.ropensci.org/targets/)
- for different cases (data, settings, etc.), we need to be able to:
   - fit a basic `corHMM` model
   - possibly: compute Wald and/or likelihood profile confidence intervals
   - specify priors and do MCMC computation
   - plot/interpret results

## contents

### documents

- `ag_bayesdiag.rmd`: Bayes diagnostics for MCMC fits ([link](http://www.math.mcmaster.ca/bolker/AG/ag_bayesdiag.html))
- `ag_supp.rmd`: supplementary material ([link](http://www.math.mcmaster.ca/bolker/AG/ag_supp.html))
- `ag_tech.rmd`: tech overview document (????)
- `ag_model.rmd`: mostly obsolete
- `corHMM.bib`: bibliography for paper/supp

### input files

- `contr*.csv`: definitions of contrast matrices for different cases. Use `source("R/utils.R"); tar_load(<R_object>); iplot(<R_object>)` to view them.
    - `contr_full.csv`: contrast matrix for full (24-parameter) model (R object `contrast_mat_full`)
	- `contr.csv`: 12-parameter model, including contrasts for *net* gain of AGs (gain - loss) (R object `contrast_mat_0`)
	- `contr_invertible.csv`: 12-parameter model with gain, loss, and including simple (1-to-1) contrasts for transitions in spawning mode and male care

### data files (`data/` directory)

- `accessTree.R`: extract phylogenies via `fishtree` package, given data
- `cleanTraitData.csv`: trait data in full format
- `binaryTraitData.csv`: trait data collapsed to binary (0/1) for traits
- `treeBlock.rds`: tree block (100 imputed trees)
- `treeSingle.rds`: phylogeny only for complete genetic data
- `dataReferencecs.csv`: list of references for each species
- `referenceList.rtf`: full reference list corresponding to the above csv

### other folders

- `bayestraits`: explorations of `BayesTraits`
   - `test_btw.R`: preliminary tests of `btw` package (BayesTraits R interface)
   - `BTrefs.xls`: papers using BayesTraits with some added information about what priors were used in the papers, how inference was done, etc..
- `old`: Miscellaneous old explorations etc.
- `twotraits`: code for two-trait analyses (comparisons with LE's thesis)


google drive [link](https://drive.google.com/drive/folders/1S5KwLDQavshwS8i0e9_g1jRiVUw8rnLO?usp=sharing)

## BayesTraits review

`BTrefs.xls` has all articles from 2020-2022 referencing Pagel & Meade (99 articles)

## installing BayesTraits

Linux (from source, V3)

```
git clone ://github.com/bbolker/BayesTraitsV3
cd BayesTraitsV3/src
## sudo apt install libnlopt-dev
make threaded
mv BayesTraitsv3_threaded ../../BayesTraitsV3
```

from binary:

```r
download.file("http://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/Files/BayesTraitsV4.0.0-Linux.tar.gz", destfile = "btv4.tgz")
untar("btv4.tgz")
unlink("btv4.tgz")
```

test:

```r
library("btw")
bayestraits(data = primate.discrete1,
	tree = primate.tree100,
	command = c("1", "1", "run"))
```

##



```
git clone git@github.com:bbolker/btw.git
```
