Here are just a few results from running cluster_se_sims.R on a local multicore
server. Doing the simulation method to learn about false positive rates (using
"coverage" of confidence interval to measure that here) and also power. Then
comparing the simulation approach to an analytic approach for power so that we
can explore power more quickly once we are satisfied that we are not going to
make too many errors with the analytic standard errors and appeals to the CLT.

So, using `cluster_se_sims.R` we see that the coverage is fine but that power is low to detect an effect of .25 sds.


Simulation based power and error rate here for 7 clusters with 1000 people in each and randomly assigning 1 to treatment:

```
r$> results <- d2_sims %>% group_by(estimator_label) %>% summarize(pow=mean(p.value < .05),
        coverage = mean(estimand <= conf.high & estimand >= conf.low),.groups= "drop")

r$> results
# A tibble: 1 x 3
  estimator_label   pow coverage
  <chr>           <dbl>    <dbl>
1 CR2 SEs         0.096        1
```

Analytic power here. The clusterPower package requires at least two clusters assigned to treatment. Here using the same ICC as above but much larger clusters.

```
r$> crtpwr.2mean(alpha=.05,m = 2, n = 100000, d = .25 , icc= .08, varw=1, cv = 0, power=NA)
 power
0.08051197

```
