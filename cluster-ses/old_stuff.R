
### This builds on https://declaredesign.org/blog/2018-10-16-few-clusters.html
num_clus <- 7
## Using 1000 instead of 100,000 for speed.
num_people_clus <- 10000
sd_cluster <- 1
rho <- 1
sd_i_1 <- .5
sd_i_0 <- 1
control_mean <- 0
treatment_mean <- .25
tau <- .25
set.seed(12345) ## same basic population for each combo
## thepop <- declare_population(
##     clusters = add_level(N = num_clus, e_c = rnorm(N) * sd_cluster, cluster_size = num_people_clus),
##     indiv = add_level(N = num_people_clus,
##         e_i0 = rnorm(N) * sd_i_0,
##         e_i1 = rnorm(n = N, mean = rho * scale(e_i0), sd = sqrt(1 - rho^2)) * sd_i_1)
## )
thepop <- declare_population(
    N = num_people_clus*num_clus,
    clusters = rep(1:num_clus,each=num_people_clus),
    e_0 = draw_normal_icc(mean = 0, N = N, clusters = clusters, sd = NULL, sd_between = NULL, total_sd = NULL, ICC = .1),
    e_1 = draw_normal_icc(mean = 0, N = N, clusters = clusters, sd = .75*sd(e_0), sd_between = NULL, total_sd = NULL, ICC = .1)
    )
potential_outcomes <- declare_potential_outcomes(Y ~ (1 - Z) * (control_mean + e_0) + Z * (tau*sd(e_0) + e_1))
estimand <- declare_estimand(ATE = mean(Y_Z_1 - Y_Z_0))
assignment_m1 <- declare_assignment(m=1,clusters = clusters)
reveal <- declare_reveal(Y, Z)
olsest <- declare_estimator(Y ~ Z, estimand = estimand, model = lm_robust, label="hc2(naive)")
cr2est <- declare_estimator(Y ~ Z, estimand = estimand, model = lm_robust, clusters = clusters, label="cr2")
cr0est <- declare_estimator(Y ~ Z, estimand = estimand, model = lm_robust, clusters = clusters,se_type="CR0", label="cr0")
stata_est <- declare_estimator(Y ~ Z, estimand = estimand, model = lm_robust, clusters = clusters,se_type="stata", label="stata")
## perm_test_fn <- function(data) {
## This doesn't work yet
##     obstz <- coef(lm(Y~Z,data=data))[["Z"]]
##     tzdist <- sapply(1:100,function(i){ ## small number of reps
##         data_c <- data %>% group_by(clusters) %>% summarize(Z=unique(Z)) %>% mutate(newZ=sample(Z))
##         data <- inner_join(data,data_c,by="clusters")
##         lmtz <- lm(Y~newZ,data=data)
##         tz <- coef(lmtz)[["newZ"]]
##         return(tz)
##     }
##         thep <- 2*min( mean(tzdist >= obstz),mean(tzdist <= obstz) )
##   return(data.frame(statistic = NA, p.value = thep))
## }
## permtest <- declare_test(handler=label_test(perm_test_fn), label="coin")
cluster_two_arm_design <- thepop +   assignment_m1 + potential_outcomes +   estimand +
    reveal + stata_est + cr0est + cr2est + olsest
dat1 <- draw_data(cluster_two_arm_design)
with(dat1,table(Z,clusters))
ICCbare(y=Y,x=clusters,data=dat1)

ptm <- proc.time()
plan(multicore)
sims <- simulate_design(cluster_two_arm_design,sims=100)
proc.time() - ptm
save(sims,file="sims.rda")
plan(sequential)
res <- sims %>% group_by(estimator_label) %>% summarize(fpr = mean(p.value < .05),
    coverage = mean(estimand <= conf.high & estimand >= conf.low, na.rm = TRUE),.groups="drop")
res
