#

## BMB for 2021-06-02

- get multi-chain MCMC working

## missing data

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
