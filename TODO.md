
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

- finish naming stuff
- fix missing data printing bug
- set up makefile
- tree blocks?
- gain/loss ratio?
   - have a parameter set and you've identified gain/loss pairs: R1, R2, ... Rn
   - right now we're passing a full parameter set and sticking parameterwise priors on top
   - but we can fairly easily add priors for (R1[1] - R1[2]), etc. : in limit we could say
       R1[gain] - R1[loss] ~ dnorm(0, sd=0.0001) -> enforcing symmetry (Sigal??)
	   differences in log-hazard rates are ("just") logs of proportional differences in gain vs loss
	   for example, we could specify (lower bound, middle, upper bound) [loss rate is 50% higher than gain,
	           middle = loss 10% higher than gain, upper = gain 20% higher than loss] - range +/- 2 or 3 SD

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