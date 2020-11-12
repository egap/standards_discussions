# The Repository for "Estimating Average Treatment Effects in Block-randomized experiments: A conversation at The No Rules of Thumb Cafe"

This repository (or subdirectory of a repository) contains the source for the standards discussion about choices involves in the estimation of average treatment effects in block-randomized experiments.

An [html preview is here](https://htmlpreview.github.io/?https://github.com/egap/standards_discussions/blob/jake_block_rand/block_rand/block_rand.html)

The pdf version can be downloaded directly.

## Provisional Precis:

Estimating the average treatment  effect in block-randomized trials requires the  use of weights. One set of weights produces  an  unbiased  estimator of the ATE  but another set  is more precise. And the extent  to  which the biased estimator  is  more  precise depends  on  the design and  data.

We  discuss the bias versus precision/variance/tradeoff  and show  by example how  a researcher can work with his  or her own data using simulation in DeclareDesign to  learn about this tradeoff.

We also remind analysts that "fixed effects" are the same as weights, and we show four different ways to produce the same estimator as the common fixed effects or least-squares dummy variable approach to the analysis of block-randomized trials.

We create some fake experiments to dramatize situations where the benefits of an unbiased estimator (using the block-size weights) outweigh the precision benefits of using the biased-estimator and also some situations where one would prefer the biased yet more precise estimator.

The overall message of this essay is that analysts need not rely on rules of thumb but can fairly easily create simulations that will reveal the bias-vs-variance tradeoffs for their own experimental designs.
