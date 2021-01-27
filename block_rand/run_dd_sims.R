## Run the simulations
library(here)
source(here::here("dd_sims.R"))
library(future)
library(future.apply)

plan(multicore)
set.seed(12345)
simvec <- c(1, 10000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
stopifnot(length(simvec) == length(designs[[1]]))
dd_sims <- simulate_designs(designs, sims = simvec)
plan(sequential)
save(dd_sims, file = here::here("dd_sims.rda"))
