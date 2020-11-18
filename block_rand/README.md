# The Repository for "Estimating Average Treatment Effects in Block-randomized experiments: A conversation at The No Rules of Thumb Cafe"

This repository (or subdirectory of a repository) contains the source for the standards discussion about choices involves in the estimation of average treatment effects in block-randomized experiments.

An [html preview is here](https://htmlpreview.github.io/?https://github.com/egap/standards_discussions/blob/jake_block_rand/block_rand/block_rand.html)

The pdf version can be downloaded directly.


## Provisional Precis:

Estimating the average treatment  effect in block-randomized trials requires
the  use of weights. One type of weights (weighting by size of block) produces
an  unbiased  estimator of the ATE as defined by the simple average of the
individual level differences in potential outcomes,  but another type of
weights (weighting by both size of block and proportion treated, often known as
"fixed effects" based weighting or an approach using "fixed effects") is more
precise. And the extent  to  which the biased estimator  is  more  precise
depends  on  the design and  data.

We  discuss the bias versus precision/variance/tradeoff  and show  by example how  a researcher can work with his  or her own data using simulation in DeclareDesign to  learn about this tradeoff.


We create some fake experiments to dramatize situations where the benefits of an unbiased estimator (using the block-size weights) outweigh the precision benefits of using the biased-estimator and also some situations where one would prefer the biased yet more precise estimator.

The overall message of this essay is that analysts need not rely on rules of thumb but can fairly easily create simulations that will reveal the bias-vs-variance tradeoffs for their own experimental designs.

## To Reproduce the document

The `Makefile` contains the relationships among the files. And the directory
includes a `renv.lock` file. If R is started from within this directory, the
`renv` system for package organization should install itself and all needed
packages.

**NOTE** Some of the steps are time consuming and may not be necessary in order
to get the main point of this piece --- that is, that when confronted by a new
block-randomized design, it might be useful to simulate in order to reason
about the bias-variance tradeoff involved in the choice of weights.

That said, here are the files and their purposes:

 - `reslst.rda` Contains the results of a search for experimental designs (in `findopts.R`) that make the block-size weighting scheme best or precision-weighting scheme best (or that make them equivalent).
 - `dd_sims.rda` Contains the results of simulating a few of those selected designs many times and using different estimators on those designs. This file come from `run_dd_sims.R`, which is a short wrapper around a few other files (see the Makefile)
 - `g_all.png` is the image in the essay. It arises from `summarize_dd_sims.R`

