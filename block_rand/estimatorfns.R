# Functions to calculate and directly apply block-size weights and precision weights for each simulation.

nbwt_est_fun <- function(data) {
  ## Block size weight
  data <- data %>%
    group_by(block) %>%
    mutate(
      newpib = mean(Z),
      newnbwt = (Z / newpib) + ((1 - Z) / (1 - newpib))
    )
  obj <- lm_robust(Y ~ Z, data = data, weights = newnbwt)
  res <- tidy(obj) %>% filter(term == "Z")
  return(res)
}

hbwt_est_fun1 <- function(data) {
  ## Precision weight
  data <- data %>%
    group_by(block) %>%
    mutate(
      newpib = mean(Z),
      newnbwt = (Z / newpib) + ((1 - Z) / (1 - newpib)),
      newhbwt = newnbwt * (newpib * (1 - newpib))
    )
  ## Not necessary since pib is same for all blocks and all assignments of Z here
  obj <- lm_robust(Y ~ Z, data = data, weights = newhbwt)
  res <- tidy(obj) %>% filter(term == "Z")
  return(res)
}

hbwt_est_fun2 <- function(data) {
  ## Precision weight
  ## Not necessary since pib is same for all blocks and all assignments of Z here
  data <- data %>%
    group_by(block) %>%
    mutate(
      newpib = mean(Z),
      newnbwt = (Z / newpib) + ((1 - Z) / (1 - newpib)),
      newhbwt = newnbwt * (newpib * (1 - newpib))
    )
  obj <- lm_robust(Y ~ Z, data = data, weights = newhbwt)
  res <- tidy(obj) %>% filter(term == "Z")
  return(res)
}

hbwt_cent_est_fun <- function(data) {
  ## Precision weight
  newdat <- data %>%
    group_by(block) %>%
    mutate(ycent = (Y - mean(Y)), zcent = (Z - mean(Z)))
  obj <- lm_robust(ycent ~ zcent, data = newdat)
  res <- tidy(obj) %>% filter(term == "zcent")
  return(res)
}


# Define estimators that can be repeated in the simulations
estnowtHC2 <- declare_estimator(Y ~ Z, estimand = "ATE", model = lm_robust, label = "Naive: Ignores Blocks")
estnbwt1 <- declare_estimator(Y ~ Z, estimand = "ATE", model = difference_in_means, blocks = block, label = "Block1: Diff Means bwt1")
estnbwt2 <- declare_estimator(handler = tidy_estimator(nbwt_est_fun), estimand = "ATE", label = "Block2: OLS bwt")
estnbwt3 <- declare_estimator(Y ~ Z, model = lm_robust, weights = (Z / Z_cond_prob) + ((1 - Z) / (Z_cond_prob)), estimand = "ATE", label = "Block3: OLS bwt")
esthbwt1 <- declare_estimator(Y ~ Z + block, estimand = "ATE", model = lm_robust, label = "Precis1: pwt Fixed Effects")
esthbwt2 <- declare_estimator(handler = tidy_estimator(hbwt_est_fun1), estimand = "ATE", label = "Precis2: OLS pwt")
esthbwt3 <- declare_estimator(handler = tidy_estimator(hbwt_est_fun2), estimand = "ATE", label = "Precis3: OLS pwt")
esthbwt4 <- declare_estimator(handler = tidy_estimator(hbwt_cent_est_fun), estimand = "ATE", label = "Precis4: OLS Center")
esthbwt5 <- declare_estimator(Y ~ Z, fixed_effects = ~block, estimand = "ATE", model = lm_robust, label = "Precis5: pwt Fixed Effects")

## Control diagnosands a bit
my_diagnosands <- declare_diagnosands(
  bias = mean(estimate - estimand),
  rmse = sqrt(mean((estimate - estimand)^2)),
  absdev = mean(abs(estimate - estimand)),
  power = mean(p.value < alpha),
  coverage = mean(estimand <= conf.high & estimand >= conf.low),
  mean_estimate = mean(estimate),
  sd_estimate = sd(estimate),
  mean_se = mean(std.error),
  mean_estimand = mean(estimand)
)
