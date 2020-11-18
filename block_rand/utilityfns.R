## Some utility functions to make the other files easier to read.

source("estimatorfns.R")

#' Make a declare design object with certain parameters
#'
#' This function takes in some parameters and produces a declare design object.
#' @param B integer number of blocks
#' @param n_b A vector of length B with the sizes of the blocks (or length 1 if all blocks are the same size)
#' @param p_b A  vector of length B with the proportions assigned to one of two experimental arms with the blocks (or length 1 if all blocks have the same proportion assigned to treatment).
#' @param const_effects_by_block Is TRUE (default) if the potential outcome to treatment is only a function of potential outcome to control or FALSE if there is some extra noise added to the potential outcome to treatment (it is no uncommon for treatment to change variance as well as location).

designer2 <- function(B, n_b, p_b, tau_b, const_effects_by_block = TRUE) {
  if (length(n_b) == 1) {
    n_b <- rep(n_b, B)
  }
  if (length(tau_b) == 1) {
    tau_b <- rep(tau_b, B)
  }
  if (length(p_b) == 1) {
    p_b <- rep(p_b, B)
  }
  stopifnot(length(p_b) == B)
  stopifnot(length(tau_b) == B)
  stopifnot(length(n_b) == B)
  set.seed(12345) ## same basic population for each combination of parameters
  U <- declare_population(
    block = add_level(
      N = B,
      prob = p_b,
      tau = tau_b
    ),
    indiv = add_level(N = n_b, e = rnorm(N))
  )
  if (!const_effects_by_block) {
    Y <- declare_potential_outcomes(
      Y_Z_0 = e,
      Y_Z_1 = e + tau + rnorm(N, mean = 0, sd = 1)
    )
  } else {
    Y <- declare_potential_outcomes(
      Y_Z_0 = e,
      Y_Z_1 = e + tau
    )
  }
  Q <- declare_estimand(ATE = mean(Y_Z_1 - Y_Z_0))
  Z <- declare_assignment(blocks = block, block_prob = p_b) # c(.5, .7, .9))
  R <- declare_reveal(Y, Z)
  ddesign <- U + Y + Z + R  +  Q + estnbwt1 + estnbwt2 + estnbwt3 + esthbwt1 + esthbwt2 + esthbwt3 + esthbwt4 + esthbwt5 + estnowtHC2
  return(ddesign)
}

dat_from_design <- function(B, n_b, p_b, tau_b, const_effects_by_block = TRUE) {
  thedes <- designer2(B = B, n_b = n_b, p_b = p_b, tau_b = tau_b, const_effects_by_block = const_effects_by_block)
  set.seed(123456) ## same data each time --- not used in declaredesign when we want to vary Z.
  thedat <- draw_data(thedes)
  ## with(dddat1, table(block, Z))
  return(thedat)
}

#' Summarize the design by block
design_summary <- function(dat) {
  within_block <-
    dat %>%
    group_by(block) %>%
    summarise(
      ATE_b = mean(Y_Z_1 - Y_Z_0),
      ATE_b_sd = sd(Y_Z_1 - Y_Z_0),
      ATE_b_est = mean(Y[Z == 1]) - mean(Y[Z == 0]),
      n_b = n(), m_b = sum(Z), p_b = mean(Z),
      ## GG 3.4
      ate_b_var = (1 / n_b) * (((n_b - m_b) * var(Y_Z_1) / m_b) + (m_b * var(Y_Z_0) / (n_b - m_b)) + 2 * cov(Y_Z_1, Y_Z_0)),
      ate_b_est_var =  (var(Y[Z == 1]) / (m_b))  + (var(Y[Z == 0]) / (n_b - m_b)),
      h_b = m_b * (n_b - m_b)/n_b,
      sample_wt = n_b,
      fe_wt1 = p_b * (1 - p_b) * n_b,
      # fe_wt2 = (2 * (m_b * (n_b - m_b))/n_b), ## same as fe_wt
      inv_var = 1 / ate_b_var,
      inv_var_est = 1 / ate_b_est_var,
      varY = var(Y),
      varYt = var(Y[Z==1]),
      varYt = var(Y[Z==0]),.groups='drop'
    ) %>%
    # divide by the sum of the weights
    mutate(
      sample_wt = sample_wt / sum(sample_wt),
      fe_wt = fe_wt1 / sum(fe_wt1),
      ## fe_wt2 = fe_wt2/ sum(fe_wt2),
      inv_var_wt = inv_var / sum(inv_var),
      inv_var_est_wt = inv_var_est / sum(inv_var_est)
    ) %>% ungroup()
  return(within_block)
}

overall_summary <- function(sum_dat) {
  res <- sum_dat %>% summarize(
    var_swt = sum(ate_b_var * sample_wt^2),
    var_fwt = sum(ate_b_var * fe_wt^2),
    var_iwt = sum(ate_b_var * inv_var_wt^2),
    estvar_swt = sum(ate_b_est_var * sample_wt^2),
    estvar_fwt = sum(ate_b_est_var * fe_wt^2),
    estvar_iwt = sum(ate_b_est_var * inv_var_est_wt^2),
    ate = sum(ATE_b * sample_wt),
    est_swt = sum(ATE_b_est * sample_wt),
    est_fwt = sum(ATE_b_est * fe_wt),
    est_iwt = sum(ATE_b_est * inv_var_est_wt),
    bias_swt = ate - est_swt,
    bias_fwt = ate - est_fwt,
    bias_iwt = ate - est_iwt,
    mse_swt = bias_swt^2 + estvar_swt,
    mse_fwt = bias_fwt^2 + estvar_fwt,
    mse_iwt = bias_iwt^2 + estvar_iwt
  )
  return(res)
}

lin_model_versions <- function(data){
    fe1_fit <- lm_robust(Y~Z,fixed_effects = ~block,data=data)
    fe2_fit <- lm_robust(Y~Z+block,data=data)
    bs_fit <- difference_in_means(Y~Z,blocks=block,data=data)
    naive_fit <- lm_robust(Y~Z,data=data)
    res <- data.frame(est_bs =   coef(bs_fit)[["Z"]],
        est_fe1 =  coef(fe1_fit)[["Z"]],
        est_fe2 =   coef(fe2_fit)[["Z"]],
        est_naive = coef(naive_fit)[["Z"]],
        estvar_bs = vcov(bs_fit)["Z","Z"],
        estvar_fe1 = vcov(fe1_fit)["Z","Z"],
        estvar_fe2 = vcov(fe2_fit)["Z","Z"],
        estvar_naive = vcov(naive_fit)["Z","Z"])
    res$ate <- with(data,mean(Y_Z_1 - Y_Z_0))
    res  <- res %>% mutate(bias_bs =  ate - coef(bs_fit)[["Z"]],
        bias_fe2 = res$ate - coef(fe1_fit)[["Z"]],
        bias_fe1 = res$ate - coef(fe2_fit)[["Z"]],
        mse_bs =  bias_bs^2 + estvar_bs,
        mse_fe1 = bias_fe1^2 + estvar_fe1,
        mse_fe2 = bias_fe2^2 + estvar_fe2)
    return(res)
}

design_to_summary <- function(B, n_b, p_b, tau_b, const_effects_by_block = TRUE,reps=1) {
    ## Notice that mse means something difference when reps==1 versus reps>1 (basically reps==1 is a rough estimate of the mse)
    thedat <- dat_from_design(B = B, n_b = n_b, p_b = p_b, tau_b = tau_b, const_effects_by_block = const_effects_by_block)
    ## Repeat experimental assignment reps times.
    if(reps==1){
        sum_dat <- overall_summary(design_summary(thedat)) %>% dplyr::select(matches("ate|fwt|swt"))
        sum2 <- lin_model_versions(thedat)
        sum3 <- cbind(sum_dat,sum2)
        return(sum3 %>% relocate(matches("mse")))
    } else {
        repsreslst <- lapply(1:reps,function(i){
            thedat <- thedat %>% group_by(block) %>% mutate(Z=sample(Z)) %>% ungroup()
            sum_dat <- overall_summary(design_summary(thedat)) %>% dplyr::select(matches("ate|fwt|swt"))
            sum2 <- lin_model_versions(thedat)
            sum3 <- cbind(sum_dat,sum2)
            return(sum3) })
        repsresdt <- do.call("rbind",repsreslst)
        names(repsresdt) <- make.names(names(repsresdt),unique=TRUE)
        stopifnot(all.equal(repsresdt$ate,repsresdt$ate.1))
        sum4 <- repsresdt %>% summarize_all(.,mean)
        ## Notice that this MSE is an average of divergencies from the truth
        sum4$mse_swt <- with(repsresdt, mean((ate - est_swt)^2 ))
        sum4$mse_fe1 <- with(repsresdt, mean((ate - est_fe1)^2 ))
        return(sum4 %>% relocate(matches("mse")))
    }
}


## Verify calcs, check for errors in liue of better testing
dat1 <- dat_from_design(B = 3, n_b = 100, p_b = c(.5, .7, .9), tau_b = c(4, 2, 0))
dat1_b <- design_summary(dat1)
dat1_overall <- overall_summary(dat1_b)
dat1_b %>% dplyr::select(block,n_b,m_b,ate_b_est_var,sample_wt,fe_wt1,fe_wt,h_b)
dat1_sum <- design_to_summary(B = 3, n_b = 100, p_b = c(.5, .7, .9), tau_b = c(4, 2, 0))
## Verify that the block-size weights are correct
e1_dat1_b1 <- difference_in_means(Y~Z,data=dat1,subset=block==1)
stopifnot(all.equal(vcov(e1_dat1_b1)[1,1], filter(dat1_b,block==1) %>% select(ate_b_est_var) %>% as.numeric() ))
e1_dat1 <- difference_in_means(Y~Z,blocks=block,data=dat1)
stopifnot(all.equal(vcov(e1_dat1)[1,1], dat1_overall$estvar_swt))
blkvar1 <- block_estimator(Y=dat1$Y,Z=dat1$Z,B=dat1$block)
stopifnot(all.equal(blkvar1$var_est,dat1_overall$estvar_swt))

dat2_sum <- design_to_summary(B = 3, n_b = c(50,150,100), p_b = c(.5, .7, .9), tau_b = c(4, 2, 0))
dat3_sum <- design_to_summary(B = 3, n_b = c(50,150,100), p_b = c(.5, .7, .9), tau_b = c(4, 2, 0), reps=10)
dat2_sum$mse_swt
dat2_sum$mse_fe1
dat3_sum$mse_swt
dat3_sum$mse_fe1
