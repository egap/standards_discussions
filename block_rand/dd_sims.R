## Simulate the designs and estimators to learn more about their properties
## This file designed to be called by run_dd_sims.R
## It uses a series of parameters for designs stored in the reslst.rda file.
## reslst arises from findopts.R (as detailed in the Makefile)

library(here)
library(knitr)
library(DeclareDesign)
library(tidyverse)
library(blkvar)
library(bfe)

source(here::here("utilityfns.R"))

source(here::here("estimatorfns.R"))

load(file = here::here("reslst.rda"))

## The columns of each element of reslst refer to a specific argument to the designer2 function such as proportion assigned to treatment in a block, etc.
lapply(reslst, round, 2)
reslstb <- lapply(reslst, as.matrix)

dd_parms0 <- tibble(
  p_b = list(
    reslstb[[1]]["hjn", 1:3],
    reslstb[[1]]["hjn1", 1:3],
    reslstb[[2]]["hjn", 1:3],
    reslstb[[2]]["hjn1", 1:3],
    reslstb[[3]]["hjn", 1:3],
    reslstb[[3]]["hjn1", 1:3],
    reslstb[[4]]["hjn", 1:3],
    reslstb[[4]]["hjn1", 1:3],
    reslstb[[5]]["hjn", 1:3],
    reslstb[[6]]["hjn", 1:3]
  ),
  tau_b = list(
    reslstb[[1]]["hjn", 4:6],
    reslstb[[1]]["hjn1", 4:6],
    reslstb[[2]]["hjn", 4:6],
    reslstb[[2]]["hjn1", 4:6],
    reslstb[[3]]["hjn", 4:6],
    reslstb[[3]]["hjn1", 4:6],
    reslstb[[4]]["hjn", 4:6],
    reslstb[[4]]["hjn1", 4:6],
    reslstb[[5]]["hjn", 4:6],
    reslstb[[6]]["hjn", 4:6]
  ),
  n_b = list(
    trunc(reslstb[[1]]["hjn", 7:9]),
    trunc(reslstb[[1]]["hjn1", 7:9]),
    trunc(reslstb[[2]]["hjn", 7:9]),
    trunc(reslstb[[2]]["hjn1", 7:9]),
    trunc(reslstb[[3]]["hjn", 7:9]),
    trunc(reslstb[[3]]["hjn1", 7:9]),
    trunc(reslstb[[4]]["hjn", 7:9]),
    trunc(reslstb[[4]]["hjn1", 7:9]),
    trunc(reslstb[[5]]["hjn", 7:9]),
    trunc(reslstb[[6]]["hjn", 7:9])
  )
)

dd_parms1 <- expand_grid(
  p_b = list(c(.5, .5, .5), c(.5, .7, .9), c(.05, .30, .40)),
  tau_b = list(c(0, 0, 0), c(.1, .1, .1), c(4, 2, 0), c(.25, .2, 0)),
  n_b = list(c(100, 100, 100), c(250, 20, 60), c(140, 80, 20))
)

dd_parms2 <- rbind(dd_parms0, dd_parms1)
dd_parms2$B <- 3
dd_parms <- rbind(dd_parms2, dd_parms2)
dd_parms$const_effects_by_block <- rep(c(TRUE, FALSE), each = nrow(dd_parms2))

save(dd_parms, file = "dd_parms.rda")

# Expand design and simulate
designs <- expand_design(
  designer = designer2,
  B = dd_parms[["B"]],
  n_b = dd_parms[["n_b"]],
  p_b = dd_parms[["p_b"]],
  tau_b = dd_parms[["tau_b"]],
  const_effects_by_block = dd_parms[["const_effects_by_block"]],
  expand = FALSE
)
