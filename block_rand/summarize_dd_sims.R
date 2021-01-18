## Visualizing the results of four designs using either block-size or precision/fixed-effect weighting

library(tidyverse)
library(DeclareDesign)
library(data.table)
library(here)
library(ggthemes)
library(ggExtra)
library(gridExtra)

load(here("dd_sims.rda"))
load(here("dd_parms.rda"))
dd_simsdt <- data.table(dd_sims)
dd_parmsdt <- data.table(dd_parms)
dd_parmsdt[, design_label := paste("design_", 1:nrow(dd_parms), sep = "")]
dd_parmsdt[, p_b_2 := lapply(p_b, round, 2)]
dd_parmsdt[, tau_b_2 := lapply(tau_b, round, 2)]
setkey(dd_simsdt, design_label)
setkey(dd_parmsdt, design_label)
dd_sims2 <- dd_simsdt[dd_parmsdt, ]
stopifnot(nrow(dd_sims2) == nrow(dd_sims))
## Add another test of the merge
rm(dd_sims) ## free memory
dd_sims2[, p_b_lab := as.character(p_b_2)]
dd_sims2[, n_b_lab := as.character(n_b)]
dd_sims2[, tau_b_lab := as.character(tau_b_2)]

## dd_sims2[sample(.N,3),.(p_b,p_b_2,p_b_lab,n_b,n_b_lab,tau_b,tau_b_2,tau_b_lab)]

## Simplify labelling
dd_sims2[, p_b_lab2 := gsub("p[0-9] = ", "", p_b_lab)]
dd_sims2[, p_b_lab2 := gsub("c(", "(", p_b_lab2, fixed = TRUE)]
dd_sims2[, p_b_lab2 := gsub("0.", ".", p_b_lab2, fixed = TRUE)]
dd_sims2[, p_b_lab2 := gsub(", ", ",", p_b_lab2, fixed = TRUE)]

dd_sims2[, n_b_lab2 := gsub("p[0-9] = ", "", n_b_lab)]
dd_sims2[, n_b_lab2 := gsub("c(", "(", n_b_lab2, fixed = TRUE)]
dd_sims2[, n_b_lab2 := gsub("0.", ".", n_b_lab2, fixed = TRUE)]
dd_sims2[, n_b_lab2 := gsub(", ", ",", n_b_lab2, fixed = TRUE)]

dd_sims2[, tau_b_lab2 := gsub("p[0-9] = ", "", tau_b_lab)]
dd_sims2[, tau_b_lab2 := gsub("c(", "(", tau_b_lab2, fixed = TRUE)]
dd_sims2[, tau_b_lab2 := gsub("0.", ".", tau_b_lab2, fixed = TRUE)]
dd_sims2[, tau_b_lab2 := gsub(", ", ",", tau_b_lab2, fixed = TRUE)]

dd_sims2[, ce_lab := as.character(const_effects_by_block)]

## dd_sims2[sample(.N,3),.(p_b,p_b_2,p_b_lab2,n_b,n_b_lab2,tau_b,tau_b_2,tau_b_lab2,ce_lab)]

new_res <- dd_sims2[, .(
  powerp = mean(p.value <= .05, na.rm = TRUE),
  powerci = mean(0 > conf.high | 0 < conf.low, na.rm = TRUE),
  meanp = mean(p.value, na.rm = TRUE),
  coverage = mean(estimand <= conf.high & estimand >= conf.low, na.rm = TRUE),
  mean_estimate = mean(estimate, na.rm = TRUE),
  bias = mean(estimate - estimand, na.rm = TRUE),
  mse = mean((estimate - estimand)^2, na.rm = TRUE),
  themad = mean(abs(estimate - estimand), na.rm = TRUE),
  sd_estimate = sd(estimate, na.rm = TRUE),
  themean_se = mean(std.error, na.rm = TRUE),
  mean_estimand = mean(estimand, na.rm = TRUE),
  nestimand = length(unique(estimand)),
  nsims = .N
), by = list(n_b_lab2, p_b_lab2, tau_b_lab2, ce_lab, estimator_label)]

new_res[, thepower := pmax(powerp, powerci)]
new_res[, mse2 := bias^2 + sd_estimate^2]

new_res[, wttype := ifelse(grepl("^Block", estimator_label), "Block-size",
  ifelse(grepl("Naive", estimator_label), "Naive", "Precision")
)]
new_res[, wttype := fct_relevel(wttype, "Block-size", "Precision", "Naive")]

tab1 <- with(new_res, table(estimator_label, wttype, exclude = c()))
stopifnot(max(unique(tab1)) == nrow(dd_parms))

## Look at the different parameters
table(new_res$n_b_lab2, exclude = c())
table(new_res$p_b_lab2, exclude = c())
table(new_res$tau_b_lab2, exclude = c())


## Create a data set that records when one or another of the weighting approaches is best
comparison1 <- new_res %>%
  filter(wttype != "Naive") %>%
  group_by(n_b_lab2, p_b_lab2, tau_b_lab2, ce_lab, wttype) %>%
  summarize(minmse = min(mse), minmad = min(themad), minbias = min(bias), .groups = "drop")
comparison1 <- comparison1 %>%
  group_by(n_b_lab2, p_b_lab2, tau_b_lab2, ce_lab) %>%
  mutate(allequal = length(unique(round(minmse, 10))) == 1) %>%
  ungroup()
comparison1 <- comparison1 %>%
  group_by(n_b_lab2, p_b_lab2, tau_b_lab2, ce_lab) %>%
  mutate(
    bs_better = minmse[wttype == "Block-size"] < minmse[wttype == "Precision"],
    pr_better = minmse[wttype == "Block-size"] > minmse[wttype == "Precision"],
    msediff = minmse[wttype == "Block-size"] - minmse[wttype == "Precision"]
  ) %>%
  ungroup()

## Now grab the ones where the block size approach is better
bs_better <- comparison1 %>%
  filter(wttype == "Block-size" & bs_better & !allequal) %>%
  arrange(desc(abs(msediff)))
## Now grab the ones where the precision based approach is better
pr_better <- comparison1 %>%
  filter(wttype == "Precision" & pr_better & !allequal) %>%
  arrange(desc(abs(msediff)))


### Plots
dd_sims2[, new_est_lab := estimator_label]
dd_sims2[new_est_lab == "Block1: Diff Means bwt1", new_est_lab := "Block-size Weights"]
dd_sims2[new_est_lab == "Precis1: pwt Fixed Effects", new_est_lab := "Precision Weights"]

#sims_pr_better <- dd_sims2[n_b_lab2 == "(138,88,18)" & p_b_lab2 == "(.05,.34,.38)" & tau_b_lab2 == "(4,3.98,4)" & ce_lab == "TRUE" &
# estimator_label %in% c("Block1: Diff Means bwt1", "Precis1: pwt Fixed Effects"), ]

sims_pr_better <- dd_sims2[n_b_lab2 == "(140,80,20)" & p_b_lab2 == "(.05,.3,.4)" & tau_b_lab2 == "(.25,.2,0)" & ce_lab == "TRUE" &
 estimator_label %in% c("Block1: Diff Means bwt1", "Precis1: pwt Fixed Effects"), ]

sims_bs_better <- dd_sims2[n_b_lab2 == "(100,100,100)" & p_b_lab2 == "(.5,.7,.9)" & tau_b_lab2 == "(4,2,0)" & ce_lab == "TRUE" &
  estimator_label %in% c("Block1: Diff Means bwt1", "Precis1: pwt Fixed Effects"), ]

sims_same <- dd_sims2[n_b_lab2 == "(250,20,60)" & p_b_lab2 == "(.5,.5,.5)" & tau_b_lab2 == "(4,2,0)" & ce_lab == "TRUE" &
  estimator_label %in% c("Block1: Diff Means bwt1", "Precis1: pwt Fixed Effects"), ]

## sims_close <- dd_sims2[n_b_lab2 == "(149,69,74)" & p_b_lab2 == "(.58,.33,.54)" & tau_b_lab2 == "(.23,.67,.19)" & ce_lab == "FALSE" &
##   estimator_label %in% c("Block1: Diff Means bwt1", "Precis1: pwt Fixed Effects"), ]

sims_close <- dd_sims2[n_b_lab2 == "(256,53,16)" & p_b_lab2 == "(.52,.53,.53)" & tau_b_lab2 == "(.38,.79,1.96)" & ce_lab == "FALSE" &
   estimator_label %in% c("Block1: Diff Means bwt1", "Precis1: pwt Fixed Effects"), ]

g_bs_better_title <- expression(paste("Exp 1: ", n[b] == "(100,100,100), ", p[b] == "(.5,.7,.9), ", tau[b] == "(4,2,0)"))
gs_bs_better_true <- sims_bs_better[, .(ate = unique(estimand))]
g_bs_better <- ggplot(sims_bs_better, aes(x = estimate, color = new_est_lab, group = new_est_lab)) +
  geom_density(size = 2) +
  geom_point(data = gs_bs_better_true, aes(x = ate, y = 0), shape = 17, size = 3, inherit.aes = FALSE) +
  geom_text(
    data = gs_bs_better_true, aes(x = ate, y = 0, label = "True ATE"), inherit.aes = FALSE,
    nudge_y = -.1
  ) +
  theme_bw() +
  labs(
    y = "Density of Estimated ATE",
    x = "Estimated Average Treatment Effect",
    title = g_bs_better_title,
    size = 2
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  theme(legend.position = c(.2, .8), legend.background = element_blank(), legend.text = element_text(size = 10)) +
  labs(color = "Estimator") +
  scale_color_colorblind() +
  guides(color = guide_legend(override.aes = list(size = 2, linetype = 1, shape = 16)))
g_bs_better
ggsave(plot = g_bs_better, filename = "g_bs_better.pdf")


g_pr_better_title <- expression(paste("Exp 2: ", n[b] == "(140,80,20), ", p[b] == "(.05,.3,.4), ", tau[b] == "(.25,.2,0)"))
gs_pr_better_true <- sims_pr_better[, .(ate = unique(estimand))]
g_pr_better <- ggplot(sims_pr_better, aes(x = estimate, color = new_est_lab, group = new_est_lab)) +
  geom_density(size = 2) +
  geom_point(data = gs_pr_better_true, aes(x = ate, y = 0), shape = 17, size = 3, inherit.aes = FALSE) +
  geom_text(
    data = gs_pr_better_true, aes(x = ate, y = 0, label = "True ATE"), inherit.aes = FALSE,
    nudge_y = -.1
  ) +
  theme_bw() +
  labs(
    y = "Density of Estimated ATE",
    x = "Estimated Average Treatment Effect",
    title = g_pr_better_title,
    size = 2
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  theme(legend.position = c(.2, .8), legend.background = element_blank(), legend.text = element_text(size = 10)) +
  labs(color = "Estimator") +
  scale_color_colorblind() +
  guides(color = guide_legend(override.aes = list(size = 2, linetype = 1, shape = 16)))
g_pr_better
ggsave(plot = g_pr_better, filename = "g_pr_better.pdf")


g_close_title <- expression(paste("Exp 3: ", n[b] == "(256,53,16), ", p[b] == "(.52,.53,.53), ", tau[b] == "(.38,.79,1.96)"))
gs_close_true <- sims_close[, .(ate = mean(estimand))]
g_close <- ggplot(sims_close, aes(x = estimate, color = new_est_lab, group = new_est_lab)) +
  #geom_density(aes(size = (3-as.numeric(as.factor(new_est_lab)))/2),alpha=.01) +
  geom_density(size=2,alpha=.1) +
  geom_point(data = gs_close_true, aes(x = ate, y = 0), shape = 17, size = 3, inherit.aes = FALSE) +
  geom_text(
    data = gs_close_true, aes(x = ate, y = 0, label = "True ATE"), inherit.aes = FALSE,
    nudge_y = -.1
  ) +
  theme_bw() +
  labs(
    y = "Density of Estimated ATE",
    x = "Estimated Average Treatment Effect",
    title = g_close_title,
    size = 2
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  theme(legend.position = c(.2, .8), legend.background = element_blank(), legend.text = element_text(size = 10)) +
  labs(color = "Estimator") +
  #scale_color_discrete(type=grDevices::adjustcolor(colorblind_pal()(3)[2:3],alpha.f=.5))+
  #discrete_scale("color",scale_name="Estimator",palette=colorblind_pal())+
  scale_color_colorblind() +
  guides(color = guide_legend(override.aes = list(size = 2, linetype = 1, shape = 16)))
g_close
ggsave(plot = g_close, filename = "g_close.pdf")


g_same_title <- expression(paste("Exp 4: ", n[b] == "(250,20,60), ", p[b] == "(.5,.5,.5), ", tau[b] == "(4,2,0)"))
gs_same_true <- sims_same[, .(ate = mean(estimand))]
g_same <- ggplot(sims_same, aes(x = estimate, color = new_est_lab, group = new_est_lab, alpha = .5)) +
  geom_density(size = 2, alpha = .5) +
  geom_point(data = gs_same_true, aes(x = ate, y = 0), shape = 17, size = 3, inherit.aes = FALSE) +
  geom_text(
    data = gs_same_true, aes(x = ate, y = 0, label = "True ATE"), inherit.aes = FALSE,
    nudge_y = -.1
  ) +
  theme_bw() +
  labs(
    y = "Density of Estimated ATE",
    x = "Estimated Average Treatment Effect",
    title = g_same_title,
    size = 2
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9)) +
  theme(legend.position = c(.2, .8), legend.background = element_blank(), legend.text = element_text(size = 10)) +
  labs(color = "Estimator") +
  scale_color_colorblind() +
  guides(color = guide_legend(override.aes = list(size = 2, linetype = 1, shape = 16)))
g_same
ggsave(plot = g_same, filename = "g_same.pdf")

g_all <- grid.arrange(g_bs_better, g_pr_better, g_close, g_same, nrow = 2, ncol = 2)

ggsave(plot = g_all, filename = "g_all.pdf", width = 10.5, height = 10.5)
ggsave(plot = g_all, filename = "g_all.png", width = 10.5, height = 10.5)
