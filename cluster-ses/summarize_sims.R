
library(DeclareDesign)
library(tidyverse)
library(here)

load(here("cluster-ses","sims.rda"))

summary(sims$estimand)

summary_sims <- sims %>% group_by(estimator_label) %>% summarize(coverage = mean(estimand <= conf.high & estimand >= conf.low),
    pow = mean(p.value < .05),
    bias = mean( estimand - estimate),sd_estimate = sd(estimate),estimand = mean(estimand), mean_se = mean(std.error))

summary_sims
