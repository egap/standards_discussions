## What if we had 7 days and we randomly assigned 1 to treatment.

library(DeclareDesign)
library(DesignLibrary)
library(dplyr)
library(ICC)
library(coin)
library(future)
library(future.apply)
## Take a short cut to get started
d1 <- block_cluster_two_arm_designer(N_blocks = 1,
    N_clusters_in_block = 7,
    N_i_in_cluster = 1000,
    sd_block = 0,
    sd_cluster = .2,
    ate = .25)
est1 <- declare_estimator(Y~Z,clusters=clusters,model=lm_robust,try_cholesky=TRUE,label="CR2 SEs")
assign1 <- declare_assignment(clusters=clusters)
d2 <- replace_step(d1,"estimator",est1)
d3 <- replace_step(d2,"assignment",assign1)

test_dat <- draw_data(d3)

library(ICC)
ICCbare(x=clusters,y=Y,data=test_dat)
plan(multicore,workers=12)
d2_sims <- simulate_design(d2,sims=c(1,1,1000,1,1,1))
plan(sequential)

results <- d2_sims %>% group_by(estimator_label) %>% summarize(pow=mean(p.value < .05),
    coverage = mean(estimand <= conf.high & estimand >= conf.low),.groups="drop")

results

summary(d2_sims$p.value)
summary(d2_sims$estimand)

save(d2_sims, file="d2_sims.rda")

## Try analytic approach for power.

library(clusterPower)
crtpwr.2mean(alpha=.05,m = 2, n = 100000, d = .25 , icc= .08, varw=1, cv = 0, power=NA)
