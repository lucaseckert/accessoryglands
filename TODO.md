## new (5 August 2024, aim for 9 August 2024)

### important priority

* "fix prior weirdness"
    * partly done (i.e. prior scaling, but LE will check and write up)
    * fix scales (i.e. fix contrast calculation, use log vs non-log scale consistently and correctly)
    * re-run and compare results with ours
    * PLOT results: (gain contrasts, loss contrasts) x 
	                (with, w/o data [i.e. prior predictive distributions]) x
	                (rj and regular) x 
					(our priors, default BayesTraits priors)
  ['with data' + 'regular' + 'our priors' should agree with our results]
  ['w/o data' + 'regular' + 'our priors' should agree with our computed priors for contrasts]
  
### lower priority

* implement multi-chain BT runs (set seeds differently), calculate R-hat etc.
* workflow cleanup (auto-run BT)
* weight contrasts by occupancy?

## new (4 Mar 2024)

* BMB: run prior-only model to extract log-mean and log-sd for rate-level priors (also look at correlations)

don't need to

```r
prior.mean <- (lb+ub)/2
prior.sd <- (ub-lb)/(2*range)
lb = log(1),
ub = log(10 * ape::Ntip(ag_compdata$phy))
## sum(branch lengths) scaled to 1
```


* BMB: try to figure out if R-hat or any other standard MCMC diagnostics are available for BT (other than acceptance rate)
* soften description of BT vs our stuff
* LE ??maybe??: run RJMCMC
* LE: try to run with 'priors only' (RJMCMC)
* BMB: weighting contrasts
* LE: what would an example look like of weighted contrasts in a simple case where we had a smaller trait space? 

## new (18 July 2023)

* LE: figure out how to run multistate model with disallowed transitions
* BMB: fight with `btw`, maybe?

## older

## BMB

* segregate core vs pix packages
* fuss with transition plots -- label on legend?
* subscripts in output table?

## BMB (27 July 2022)

https://docs.google.com/document/d/1caLrGnvanUawtSbiyjJ8PkL6JkecZX5f/edit?usp=sharing&ouid=115557938223853461511&rtpof=true&sd=true

### Main paper

- check with Beaulieu et al. about preferred ref
- Lit review/sample of how BayesTraits is used?
   - inference?
   - default priors/mention of priors?
- caveats para (hard-coded root; 12 vs 24 parameters)? Comment on (lack of) sensitivity?


### Supplement

- PDF/docx output? (code folding only in HTML ...)
- Work harder on getting reliable CIs for MLE fits
- Roll bayesdiag into supplement (or, keep as a separate doc?)
- tech supp
   - static picture of `targets` flow?
   
### Low priority

- rerun with prior on root rather than hard value?
- Compute Bayes Factors???
    - Stepping-stone algorithm??? see Baele and Vansteelandt 2013
		
Baele, Guy, Philippe Lemey, and Stijn Vansteelandt. “Make the Most of Your Samples: Bayes Factor Estimators for High-Dimensional Models of Sequence Evolution.” BMC Bioinformatics 14, no. 1 (March 6, 2013): 85. https://doi.org/10.1186/1471-2105-14-85.


- Weight contrasts by state occupancy?
- dig out/discuss description of BayesTraits priors?
- Mention that losses are uncertain/less well estimated because AGs are absent in most of the tree?
- cleanup/code stuff
   - use `renv` for package versions?
   - priors should be specified externally to `_targets` file
   - make sure thinning params get passed through to MCMC function (but switch targets to use 10 rather than 20 to match current setup ...)
   - describe/check lower/upper bounds in `ag_model_tb` (0.1 to 100N) vs `ag_model_pscs_prior` (1 to 10N) ??
   - double-check thinning? (supposedly ran burn-in of 4000 + 80000 steps with thinning of 20; how come we have n = 8000 in the resulting chains? thinning arg didn't get passed through)
   
## LE for 2021-07-07

- try all of the following for the two trait model
  1. BayesTraits and BayesFactors - DONE
  2. Reversible Jump Models
  3. Bayesian with contrasts (with tree block)
  4. Pagels test - DONE
  5. likelihood ratio test - DONE
  6. AIC comparison - DONE
  7. BayesTraits with contrasts
- put it all together in a rmd to share with Sigal
- compute ESS for bayesian stuff - DONE
- sanity check on contrasts

## BMB for 07-07-30

- weight contrasts by state occupancy (????)
- likelihood (point) comparison, 4 way (independent, pc-only, sc-only, full {pc+sc})
- move tech/devel stuff out of `ag_model.rmd`, more text/explanation
- reorder params before Bayes pairs plots?

- consider run/comparison for spawning vs repro system?
- why is evidence apparently stronger in 2-trait model?
    - joint test on 2 parameters (gain diff, loss diff) ?
	- dilution by non-focal parameters?
- general sanity check of parameters (DONE)
- sanity check of treeblock runs (DONE)
- compute improved R-hat (DONE)
- visualize tip state distribution for non-treeblock/treeblock (DONE)

## LE for 2021-06-30

- run two trait model with care and compare to previous results
- make sure roots are set

## BMB for 2021-06-30

- think about tree blocks
- implement priors/scaling (DONE)
- shorten names? (DONE)
- parameter ordering and log-scaling for params plot (DONE)
- adjust `getdata` to allow flexibility in selected trait columns (LOW)
- dynamic branching (for primates + ag + ag w/o missing data + .. ? (LOW)

## LE for 2021-06-23

- finish contr - done
- synthesize advice on priors (gain/loss ratios)
- primate comparison mcmc
- think about why some pars are so uncertain (trait combinations)

## BMB for 2021-06-23

- shorten names?
- think about tree blocks
- adjust `getdata` to allow flexibility in selected trait columns
- dynamic branching (for primates + ag + ag w/o missing data + .. ?
- parameter ordering and log-scaling for params plot
- phylo scaling (pre- or post-hoc)

## LE for 2021-06-16

- discrete vs continuous comparison (Pagel/Meade/etc vs parameter estimation)
    - Pagel & Meade
	- Beaulieu O'Meara
	- your thesis
- list of stuff to do/intentions, classified by importance and difficulty
- figure out where the tree block data comes from ... ?
- interactions
- https://www.embarcadero.com/starthere/xe5/mobdevsetup/ios/en/installing_the_commandline_tools.html

### expanded histograms

```r
d <- as.data.frame(matrix(rnorm(6000), ncol=6))
library(tidyr)
d_long <- pivot_longer(d, cols=everything(), names_to="var")
library(ggplot2)
ggplot(d_long, aes(x=value)) + facet_wrap(~var) +
  geom_histogram(binwidth=0.1)
```

## Sigal and Jess

- timing and venue
- gain/loss priors

## BMB 

- finish naming stuff (done)
- fix missing data printing bug: partially done (no longer breaks, but also doesn't match number of states? collapse vs non-collapse?)
- keep working on `targets`:
   - dynamic 
   - target documentation?
- tree blocks?
- gain/loss ratio?
   - have a parameter set and you've identified gain/loss pairs: R1, R2, ... Rn
   - right now we're passing a full parameter set and sticking parameterwise priors on top
   - but we can fairly easily add priors for (R1[1] - R1[2]), etc. : in limit we could say
       R1[gain] - R1[loss] ~ dnorm(0, sd=0.0001) -> enforcing symmetry (Sigal??)
	   differences in log-hazard rates are ("just") logs of proportional differences in gain vs loss
	   for example, we could specify (lower bound, middle, upper bound) [loss rate is 50% higher than gain,
	           middle = loss 10% higher than gain, upper = gain 20% higher than loss] - range +/- 2 or 3 SD

* priors on branch length:

## LE for 2021-06-09

- figure out ?s for missing data
- build up contrasts by computing intermediate values, e.g.

To identify a parameter we need to know the
both the row ('from') and column ('to')
name variables as 'gain' and loss'
```
{ag0_pc0_sc0, ag0_pc0_sc1} -> "gain.sc_ag0_pc0"
{ag0_pc0_sc1, ag0_pc0_sc0} -> "loss.sc_ag0_pc0"

## difference in gain rates for sc1 vs sc0 given pc0
sc.gain.pc0 <- gain.ag_sc1_pc0 - gain.ag_sc0_pc0
sc.gain <- (sc.gain.pc1 + sc.gain.pc0)/2
sc.netgain <- sc.gain - sc.loss
```

- bivariate plot of `sc.gain` vs `sc.loss`
- histograms or violin plots of `sc.gain`, `sc.loss`, `sc.netgain`

```r
var     sample   value
sc.gain 1        ...
sc.gain 2        ...
...
sc.loss 1        ...
```

```r
ggplot(df, aes(x=var, y=value)) + geom_violin()
```

- is there a parameter contrast for combined effects of PC and SC?
- loss of power due to missing data? (either in phylog or in traits)

```r
## after introducing "?" for missing values and reading ?corHMM
corHMM:::corProcessData(data) ## then figure out what all the bits are
```

## BMB for 2021-06-02

- multi-chain MCMC, HPD regions on pairs plots (done)
- helper function for variable names?
- look into missing-value computations ...

## missing data (solved!)

missing data, mostly for parental care/mating system traits (AG info is known)

- maximum fallback: forget it all, go back to BayesTools
- or: run pairwise analyses in augmented corHMM (i.e. get confidence intervals etc.)
- or: single-trait single imputation (continuous hack via phytools or invent something. Use ancestral character estimation/mapping to MRCA of known and unknown taxa/trait combos
- or: multiple imputation? ugh
- or: see if we can roll an 'unknown state' and/or unobserved values into the machinery

## tree blocks/uncertainty

- consensus tree from tree block?
- use complete tree?

## miscellaneous/inference

- think about multimodal likelihood surface
   - MCMC algorithms
   - multiple chains?
   - run lots longer?
   - thin?
   
- gain/loss polarity

## next steps

- clear picture of ideal program
- figure out how to parameterize interactions
- compare monkey/Pagel & Meade?
- power differences?
- biologically
- weighting when averaging across states?

## machinery

- figure out a good workflow that handles slow bits? Do we need `make`, or `targets`, now?
- for now: Rmarkdown with caching, plus explicit caching
- LE: install `styler` add-in, use it to clean up formatting as you go: https://styler.r-lib.org/
- BMB: implement multi-chain/multi-start/parallel M-H (with tuning??) e.g. start with the code [here](https://cran.microsoft.com/web/packages/glmmTMB/vignettes/mcmc.html) (+ pick random starting values from the priors, use Hessian for scale; adapt scale??)  (NIMBLE package is suggested; or the Metropolis function from [BayesianTools](https://github.com/florianhartig/BayesianTools/blob/master/BayesianTools/man/Metropolis.Rd); `ramcmc` package?
- image method for rate matrices? improved parameter naming?

## Qs for corHMM authors

- what about identifiability of gain/loss?
- are they interested in the philosophical/CI-vs-hypothesis-comparison stuff?
