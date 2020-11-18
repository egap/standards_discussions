
library(here)
source(here::here("new_sims.R"))

plan(multicore)
set.seed(12345)
simvec <- c(1, 1000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
stopifnot(length(simvec) == length(designs[[1]]))
new_sims <- simulate_designs(designs, sims = simvec)
plan(sequential)
save(new_sims, file = here::here("new_sims.rda"))
