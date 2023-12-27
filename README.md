AdaptiveR
================

Functions included in `AdaptiveR` helps users implement a simple
adaptive survey design with standard Bernoulli Thompson sampling
introduced in [Lee and Green
(2023+)](https://www.dropbox.com/scl/fi/o0o62xwd278lrrtrisg9v/RCV_adaptive.pdf?dl=0&rlkey=owlw1xuqr4o4tnyx1oxl4ey3s).
It includes functions that 1) iteratively update treatment assignment
probability based on Thompson sampling algorithm as new data comes in
and, at the end of the experiment, 2) estimate mean outcomes using
Inverse Probability Weight (IPW) estimator, 3) create a coefficient plot
based on the estimated mean outcomes, and 4) create a plot visualizing
over-time posterior probability development.

### Installation

``` r
devtools::install_github('DianaDaInLee/AdaptiveR')
```

### Adaptive Experiment

To demonstrate the functionalities in `AdaptieR`, we use an example of
an adaptive survey experiment with 5 treatment arms sampling 50
observations over 10 periods. The function `run_adapt` allows users to
compute posterior probability of being best arm using standard Bernoulli
Thompson sampling. It returns matrices containing posterior
probabilities, total number of observations, total successes (outcome
equals 1), as well as survey data with treatment assignment probability
and IPW for each observation appended as new columns.

**Arguments**

- `prior_survey_data` (default = `NULL`) data.frame containing survey
  results from prior period. Must contain the following two columns:
  “arm” that contains treatment arm numbers 1 through k, and “Y” that
  contains binary survey outcome (0 or 1). Ignored in the first period
  (i.e., when `current_period` == 1).
- `prior_adapt_matrix` (default = `NULL`) result of `run_adapt()` from
  the previous batch. Ignored in the first period (i.e., when
  `current_period` == 1).
- `current_period` (default = `NULL`) numeric value indicating the
  number of current batch to be implemented (i.e., total number of
  batches previously done + 1). For example, if a user has conducted a
  survey for the first batch and is looking to update treatment
  probability for the second batch, the value should be 2.
- `arms` (default = `NULL`) Total number of arms in the survey.
- `periods` (default = `NULL`) Total number of periods to be implemented
  in the survey.
- `floor` (default = `NULL`) A numeric value between 0 and 1 indicating
  the minimum probability of treatment assignment imposed on each arm.
- `seed` (default = 1004)

#### Initial Period

In the first period, given no survey results exist, `run_adapt` will
only require `period`, `arms` and `current_period` arguments to be
filled in to generate treatment assignment probabilities as well as to
initialize matrices to store information in the future periods.

``` r
x <- run_adapt(period = 10, arms = 5, current_period = 1)
```

    ## Treatment assignment probabilities for period 1 are as follows:
    ##  0.2, 0.2, 0.2, 0.2, 0.2

``` r
names(x)
```

    ## [1] "m_posterior"       "m_trials"          "m_success"        
    ## [4] "prior_survey_data"

``` r
head(x$m_posterior, n = 3)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]  0.2  0.2  0.2  0.2  0.2
    ## [2,]   NA   NA   NA   NA   NA
    ## [3,]   NA   NA   NA   NA   NA

The function in the first period will create three
matrices—`m_posterior`, `m_trials` and `m_success`—with rows
representing number of periods and columns representing total treatment
arms. The first row (i.e., first period) in `m_posterior` will be filled
in with the treatment assignment probability while all matrices will be
entirely empty as no new data has been collected. Similarly,
`prior_survey_data` will be `NULL`.

#### Subsequent Period

In the second period, survey results as well as the `run_adapt` object
from the first period should be supplied in order to update the
treatment assignment probability. We will use a fake data
`result_batch1` in `ex_survey_data` created for demonstration purposes:

``` r
data(ex_survey_data)
result_batch1 <- ex_survey_data[ex_survey_data$period==1, c('arm', 'Y')]
head(result_batch1, n = 3)
```

    ##   arm Y
    ## 1   3 0
    ## 2   1 0
    ## 3   3 1

Re-run `run_adapt`:

``` r
x <- run_adapt(current_period = 2, prior_survey_data = result_batch1, prior_adapt_matrix = x, floor = 0.02)
```

    ## Treatment assignment probabilities for period 2 are as follows:
    ##  0.33211118464593, 0.472481800132363, 0.129199205823958, 0.0462078093977498, 0.02

``` r
head(x$m_posterior, n = 2)
```

    ##           [,1]      [,2]      [,3]       [,4] [,5]
    ## [1,] 0.2000000 0.2000000 0.2000000 0.20000000 0.20
    ## [2,] 0.3321112 0.4724818 0.1291992 0.04620781 0.02

``` r
head(x$m_success, n = 2)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    6    9    6    2    2
    ## [2,]   NA   NA   NA   NA   NA

``` r
head(x$m_trials, n = 2)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]   10   14   12    6    8
    ## [2,]   NA   NA   NA   NA   NA

``` r
head(x$prior_survey_data, n = 2)
```

    ##   arm Y period probs ipw
    ## 1   3 0      1   0.2   5
    ## 2   1 0      1   0.2   5

For each of the subsequent periods, run `run_adapt` to update treatment
assignment probabilities as well as to record new information as they
come in. Note that, to record information from the last period, run with
`current_period` set equal to the last period + 1 (i.e., 11 in the case
of 10 periods). This will allow for success, trials, and survey data
from the 10th period to be appended to the `run_adapt` object.

### Post Data Collection

Once the data collection is completed, users can use other functions to
get estimates and visualize the experiment results. Here, we work with
an example data, which is a return from `run_adapt()` after completing
all 10 periods of experiment:

``` r
data("ex_run_adapt_data")
```

#### Posterior Probability

`plot_posterior` provides ggplot-based plot of over-time posterior
probability development.

``` r
p_plot <- plot_posterior(adapt_matrix = ex_run_adapt_data)
p_plot
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### Mean Outcome

`adapt_est` estimates inverse probability weighted estimated mean
outcome using `estimatr::lm_robust` and creates a coefficient plot based
on the result.

``` r
est <- adapt_est(adapt_matrix = ex_run_adapt_data)
est$est
```

    ##       Estimate Std. Error   t value     Pr(>|t|)    CI Lower  CI Upper  DF
    ## arm1 0.6437764 0.02776316 23.188151 4.492821e-81  0.58922826 0.6983246 495
    ## arm2 0.5112921 0.11130386  4.593660 5.532374e-06  0.29260588 0.7299784 495
    ## arm3 0.4322043 0.11174801  3.867669 1.245546e-04  0.21264538 0.6517632 495
    ## arm4 0.3216712 0.12887833  2.495930 1.288708e-02  0.06845521 0.5748872 495
    ## arm5 0.2500000 0.14848930  1.683623 9.288502e-02 -0.04174702 0.5417470 495

``` r
est$est_plot
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Simulation

Users can use `simulation` function to evaluate the performance of
either static or adaptive design under settings specific to their
research. It returns several matrices that stores simulation results as
well as simulated statistics including RMSE, bias, and coverage that
helps users assess optimal parameters (number of periods, arms, sample
size). We encourage users to utilize this function to test various
settings and identify the experimental design most appropriate for their
research conditions and/or assess the relative benefit of adaptive
design to static design.[^1]

**Arguments**

- `probs` (default = `NULL`) distribution of true mean outcomes for each
  arm
- `static` (default = `FALSE`) whether the simulated experiment should
  be adaptive or static
- `periods` (default = 10) total number of batches in a single
  experiment (relevant for adaptive design only)
- `n` (default = 1000) total sample size
- `n_first` (default = `NULL`) total sample size to be allotted to the
  first batch. If `NULL`, then `n` is evenly distributed across batches
- `floor_rate` (default = 0.01) minimum sample probability to assign to
  ensure all arms receive non-zero sample in each batch
- `iter` (default = 1000) number of simulations to run

``` r
sim_out <- simulate(probs = c(0.2, 0.15, 0.1, 0.1), iter = 10)
```

    ## 
    ## -----------------------------
    ## Adaptive simulated experiment will be performed with following parameters:
    ##  Number of Arms: 4 
    ##  Number of Periods: 10 
    ##  Total Sample Size: 1000 
    ##  Sample Size for First Period: 100 
    ##  Sample Size for Remaining Periods: 100 
    ## -----------------------------
    ## 
    ##  Iteration: 1 ..2 ..3 ..4 ..5 ..6 ..7 ..8 ..9 ..10 ..

``` r
sim_out$outmat
```

    ## # A tibble: 4 × 6
    ##   term  true  best  bias   rmse  coverage
    ##   <chr> <chr> <chr> <chr>  <chr> <chr>   
    ## 1 Arm 1 0.200 1.000 -0.002 0.015 0.900   
    ## 2 Arm 2 0.150 0.000 0.006  0.054 0.800   
    ## 3 Arm 3 0.100 0.000 -0.024 0.079 0.800   
    ## 4 Arm 4 0.100 0.000 0.043  0.062 0.500

[^1]: This function is an adaptation of R scripts provided as part of
    the replication materials of [Offer-Westort, Coppock, and Green
    (2021)](https://onlinelibrary.wiley.com/doi/abs/10.1111/ajps.12597).
