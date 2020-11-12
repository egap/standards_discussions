
library(DeclareDesign)
library(tidyverse)
library(here)

load(here("cluster-ses","sims.rda"))

unique(sims$estimand)

summary_sims <- sims %>% group_by(estimator_label) %>% summarize(false_pos_error= mean( estimand <= conf.low | estimand >= conf.high),
    #coverage = 1 - false_pos_error,
    coverage = mean(estimand <= conf.high & estimand >= conf.low),
    pow = mean(p.value < .05),
    bias = mean( estimand - estimate),sd_estimate = sd(estimate),estimand = mean(estimand), mean_se = mean(std.error))

summary_sims
