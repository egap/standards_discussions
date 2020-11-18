## Expanding on block-randomized-exps.R
## Trying to figure out why intution of block vs precision weights
## does not translate to simulations

library(DeclareDesign)
library(tidyverse)
library(future)
library(future.apply)
library(blkvar)
library(bfe)
library(optimr)

source(here::here("utilityfns.R"))

fn1a <- function(x, consteffects = TRUE, optwhat = "bestfe",reps=1) {
  if (any(x[1:3] >= 1) || any(x[1:3] <= 0) || any(x[1:3] < .05) || any(x[1:3] > .95) || any(x[4:6] < 0) || any(x[4:6] > 4)) {
    res <- 9999999
    return(res)
  }
  p_b <- round(x[1:3], 2)
  n_b <- trunc(x[7:9])
  n_t <- trunc(n_b * p_b)
  n_c <- n_b - n_t
  tau_b <- round(x[4:6], 2)
  ## No more than 400 and no less than 200 units total
  if (any(n_t < 4) || any(n_c < 4) || sum(n_b) < 200 || sum(n_b) > 400) {
    res <- 9999999
    return(res)
  }
  des <- design_to_summary(B = 3, n_b = n_b, p_b = p_b, tau_b = tau_b, const_effects_by_block = consteffects,reps=reps)
  ## Maybe and probably shuffle Z and recalc the summary a couple of times.
  if (optwhat == "bestfe") {
    return(des$mse_fe1 - des$mse_swt)
  }
  if (optwhat == "bestbs") {
    return(des$mse_swt - des$mse_fe1)
  }
  if (optwhat == "same") {
    return((des$mse_swt - des$mse_fe1)^2)
  }
}

# fn1(c(.5, .7, .9, 4, 2, 0))
# fn1(c(.5, 1.7, .9, 4, 2, 0))
# fn1(c(.05, .2, .9, 4, 2, 0))
# fn1(c(.05, .2, .9, 4, 2, 0), consteffects = FALSE)
fn1a(c(.05, .2, .9, 4, 2, 0, 100, 100, 100), consteffects = TRUE)
fn1a(c(.05, .2, .9, 4, 2, 0, 100, 100, 100), consteffects = FALSE)
fn1a(c(.05, .2, .9, 4, 2, 0, 100.1, 100.2, 100.3), consteffects = FALSE,optwhat="bestfe",reps=1)
fn1a(c(.05, .2, .9, 4, 2, 0, 100.1, 100.2, 100.3), consteffects = FALSE,optwhat="bestfe",reps=10)
fn1a(c(.05, .2, .9, 4, 2, 0, 100.1, 100.2, 100.3), consteffects = FALSE,optwhat="bestfe",reps=10)

opts_fn_a <- function(fn, consteffects, optwhat) {
  ## The idea here is to do two wide stochastic searches and then use better optimizers
  theparscale <- c(1, 1, 1, 4, 4, 4, 400, 400, 400)
  ## optmethods <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "nlminb", "Rcgmin", "Rvmmin", "hjn")
  optmethods <- c("Rcgmin", "Rvmmin", "hjn")
  ## Start with a design and data scenario that is best for block-size weights
  ## opt1 <- optim(par=c(.5,.7,.9,4,2,0,333,333,334),fn=fn1a,consteffects=TRUE,optwhat="equal",method="SANN",control=list(trace=1,parscale=theparscale,temp=100,tmax=100))
  opt1 <- optim(par = c(.5, .7, .9, 4, 2, 0, 100, 100, 100), fn = fn, consteffects = consteffects,
      optwhat = optwhat, reps=8, method = "SANN", control = list(trace = 1, parscale = theparscale, temp = 100, tmax = 100))
  opt1best <- c(opt1$par, opt1$value)
  ## Add a scenario that is neutral --- different weights amount to the same thing
  opt2 <- optim(par = c(.5, .5, .5, 0, 0, 0, 100, 100, 100), method = "SANN", fn = fn1a, consteffects = consteffects,
      optwhat = optwhat, reps=8, control = list(trace = 1, parscale = theparscale, temp = 100, tmax = 100))
  opt2best <- c(opt2$par, opt2$value)
  opt4 <- opm(par = opt1$par, fn = fn1a, consteffects = consteffects, optwhat = optwhat, reps=8, method = optmethods,
      control = list(trace = 1, parscale = theparscale))
  opt4$method <- row.names(opt4)
  opt4best <- opt4[opt4$value == min(opt4$value), 1:10]
  opt5 <- opm(par = opt2$par, fn = fn, consteffects = consteffects, optwhat = optwhat, reps=8, method = optmethods,
      control = list(trace = 1, parscale = theparscale))
  opt5$method <- row.names(opt5)
  opt5best <- opt5[opt5$value == min(opt5$value), 1:10]
  return(rbind(opt1best, opt2best, opt4best, opt5best))
}


library(progressr)
parms <- expand.grid(consteffects = c(TRUE, FALSE),optwhat=c("bestfe","bestbs","same"))

plan(multicore, workers = 6)
is <- 1:6
with_progress({
  p <- progressor(along = is)
  reslst <- future_lapply(is, function(i, ...) {
    theparms <- parms[i, ]
    p(sprintf("x=%g", i))
    res <- opts_fn_a(fn = fn1a, consteffects = theparms[1], optwhat = theparms[2])
    fnm <- paste(theparms, "res_tmp.rda", collapse = "")
    save(res, file = fnm)
    return(res)
  },future.seed=12345)
})
plan(sequential)
nms <- apply(parms,1,paste0,collapse="_")
names(reslst) <- nms
save(reslst, file = here::here("reslst.rda"))




## res1 <- opts_fn(fn=fn1,consteffects=TRUE)
## res2 <- opts_fn(fn=fn1,consteffects=TRUE)
## res3 <- opts_fn(fn=fn1,consteffects=TRUE)
## res4 <- opts_fn(fn=fn1,consteffects=TRUE)
##

## parms <- data.frame(fn = c("fn1a", "fn1a", "fn2a", "fn2a"), consteffects = c(TRUE, FALSE, TRUE, FALSE))
##
## library(progressr)
##
## plan(multicore, workers = 4)
## is <- 1:4
## with_progress({
##   p <- progressor(along = is)
##   reslst <- future_lapply(is, function(i, ...) {
##     p(sprintf("x=%g", i))
##     theparms <- parms[i, ]
##     fnm <- paste(theparms[1], theparms[2], "res_tmp.rda", collapse = "")
##     res <- opts_fn(fn = eval(parse(text = theparms[1])), consteffects = theparms[2])
##     save(res, file = fnm)
##     return(res)
##   })
## })
##
## plan(sequential)
## nms <- c("pre_best_ce", "pre_best_noce", "bs_best_ce", "bs_best_noce")
## names(reslst) <- nms
## save(reslst, file = here::here("reslst.rda"))
##
##
## ## pre_best_ce <- opts_fn(fn=fn1,consteffects=TRUE)
## ## pre_best_noce <- opts_fn(fn=fn1,consteffects=FALSE)
## ## bs_best_ce <- opts_fn(fn=fn2,consteffects=TRUE)
## ## bs_best_noce <- opts_fn(fn=fn2,consteffects=FALSE)
## ## save(pre_best_ce,pre_best_noce,bs_best_ce,bs_best_noce, file=here::here("optsout.rda"))
##
##
## ## Try grid search to draw nice graphs
##
## newparms <- expand.grid(seq(.05, .95, .05), seq(0, 4, .25), seq(40, 920, 20))
