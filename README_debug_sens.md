Something is messed up with the sensitivity analysis plots in the supp.

`old_sens.png` is a screenshot of the compiled figure from a stats talk, compiled 17 Jan 2023.  This makes sense and matches with the text description in the supp figure.



## attempting to resurrect old version/git bisect

This seems to be problematic: all kinds of headaches with `targets` versions, dependencies, etc..

https://stackoverflow.com/questions/6990484/how-to-checkout-in-git-by-date:

```
git checkout `git rev-list -1 --before="Jan 15 2023" HEAD`
```

this goes back to f9e42f7

have checked this out as a `salvage` branch ...

```r
## need to downgrade targets if we don't want to rerun everything
## (back to 1.6.0 for pipeline changes, then conflict with secretbase::sha3() ?
Rscript -e 'remotes::install_version("tarchetypes", "0.7.3")'
Rscript -e 'remotes::install_version("targets", "0.14.1")'
Rscript -e 'targets::tar_make("ag_supp_html")'
```

After this, though, we get to 

> 'ag_mcmc_0' not found

I'm afraid this is going to be slow to rebuild ... instead try to stare at `git diff` to see if I can figure it out.

(But maybe it rebuilt anyway when I did `tar_make` again??)

May be able to `git bisect`, but will be very slow (rebuilding everything takes time).

In the meantime, can I figure out what changed?  General themes of what we've done:

* BayesTraits work (should be independent/irrelevant unless we touched something important by accident?)
* missing data/subsampling examples
* flow figure
* maybe contrast weighting??? (3d10c1e7)

```
git diff f9e42f7 HEAD
```
