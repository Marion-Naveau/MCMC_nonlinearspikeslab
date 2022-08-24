**n.b.:** to update this file, you can change README.Rmd then run
`make readme` in the command line.

# nonlinearspikeslab

<!-- badges: start -->
<!-- badges: end -->
<!--ts-->
* [nonlinearspikeslab](#nonlinearspikeslab)
   * [Installation](#installation)
   * [Example 1: regression on the final size parameter](#example-1-regression-on-the-final-size-parameter)
      * [View simulated data](#view-simulated-data)
      * [Fit Spike and Slab model with fixed parameters](#fit-spike-and-slab-model-with-fixed-parameters)
      * [Fit Spike and Slab model](#fit-spike-and-slab-model)
      * [Prior posterior plot](#prior-posterior-plot)
   * [Example 2: regression on the delay parameter, after reparametrisation](#example-2-regression-on-the-delay-parameter-after-reparametrisation)
      * [View simulated data](#view-simulated-data-1)

<!-- Added by: gkonkamking, at: mar. 22 févr. 2022 18:40:44 CET -->

<!--te-->

Spike and slab variable selection in nonlinear models

## Installation

Clone the repository to a local folder using a command like:

``` bash
git clone git@forgemia.inra.fr:konkam/nonlinearspikeslab.git
```

or in case the previous does not work:

``` bash
git clone https://forgemia.inra.fr/konkam/nonlinearspikeslab.git
```

Then:

1.  Open RStudio, open a new project inside the newly created folder

2.  Press ctrl + maj + B

## Example 1: regression on the final size parameter

``` r
library(nonlinearspikeslab)
library(nimble)
```

    ## nimble version 0.12.1 is loaded.
    ## For more information on NIMBLE and a User Manual,
    ## please visit https://R-nimble.org.

    ## 
    ## Attaching package: 'nimble'

    ## The following object is masked from 'package:stats':
    ## 
    ##     simulate

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.6     ✓ dplyr   1.0.8
    ## ✓ tidyr   1.2.0     ✓ stringr 1.4.0
    ## ✓ readr   2.1.2     ✓ forcats 0.5.1

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
set.seed(0)
data = SimulateDataOneVariableParam(p=10)
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

    ## Joining, by = "i"

``` r
data
```

    ## # A tibble: 300 × 28
    ##    times     i phi_i  gphi     Y    V2     V3    V4    V5     V6    V7      V8
    ##    <dbl> <int> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>
    ##  1   0       1  16.7  8.35  8.40  1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  2  22.2     1  16.7 11.8  11.8   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  3  44.4     1  16.7 14.3  14.3   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  4  66.7     1  16.7 15.6  15.6   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  5  88.9     1  16.7 16.2  16.1   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  6 111.      1  16.7 16.5  16.5   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  7 133.      1  16.7 16.6  16.7   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  8 156.      1  16.7 16.7  16.7   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  9 178.      1  16.7 16.7  16.7   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ## 10 200       1  16.7 16.7  16.7   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ## # … with 290 more rows, and 16 more variables: V9 <dbl>, V10 <dbl>, V11 <dbl>,
    ## #   phi2 <dbl>, phi3 <dbl>, alpha <dbl>, beta1 <int>, beta2 <int>, beta3 <int>,
    ## #   beta4 <int>, beta5 <int>, beta6 <int>, beta7 <int>, beta8 <int>,
    ## #   beta9 <int>, beta10 <int>

### View simulated data

``` r
PlotSimulatedData(data)
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

### Fit Spike and Slab model with fixed parameters

``` r
mcmc_out = FitOneVariableParamSpikeAndSlab_fixed_sigmagamma(data, nu_0 = 10^(-4))
```

    ## Defining model

    ## Building model

    ## Setting data and initial values

    ##   [Note] 'sigma' has initial values but is not a variable in the model and is being ignored.

    ##   [Note] 'Gamma' has initial values but is not a variable in the model and is being ignored.

    ## Running calculate on model
    ##   [Note] Any error reports that follow may simply reflect missing values in model variables.

    ## Checking model sizes and dimensions

    ## ===== Monitors =====
    ## thin = 1: alpha, mu, omega
    ## ===== Samplers =====
    ## RW sampler (2)
    ##   - omega
    ##   - mu
    ## conjugate sampler (41)
    ##   - alpha
    ##   - beta[]  (10 elements)
    ##   - phi[]  (30 elements)
    ## binary sampler (10)
    ##   - delta[]  (10 elements)

    ## Defining model

    ## Building model

    ## Setting data and initial values

    ##   [Note] 'sigma' has initial values but is not a variable in the model and is being ignored.

    ##   [Note] 'Gamma' has initial values but is not a variable in the model and is being ignored.

    ## Running calculate on model
    ##   [Note] Any error reports that follow may simply reflect missing values in model variables.

    ## Checking model sizes and dimensions

    ## Checking model calculations

    ## Compiling
    ##   [Note] This may take a minute.
    ##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.

    ## running chain 1...

    ##   [Note] 'sigma' has initial values but is not a variable in the model and is being ignored.

    ##   [Note] 'Gamma' has initial values but is not a variable in the model and is being ignored.

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|

    ## running chain 2...

    ##   [Note] 'sigma' has initial values but is not a variable in the model and is being ignored.

    ##   [Note] 'Gamma' has initial values but is not a variable in the model and is being ignored.

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|

    ## running chain 3...

    ##   [Note] 'sigma' has initial values but is not a variable in the model and is being ignored.

    ##   [Note] 'Gamma' has initial values but is not a variable in the model and is being ignored.

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|
    ## There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.

``` r
posterior_on_number_of_covariates_selected(mcmc_out, warmup_prop = 0)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
traceplot(mcmc_out)
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
posterior_predictive_estimate(data = data, mcmc_out = mcmc_out)
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

### Fit Spike and Slab model

``` r
mcmc_out = FitOneVariableParamSpikeAndSlab(data, nu_0 = 10^(-4))
```

    ## Defining model

    ## Building model

    ## Setting data and initial values

    ## Running calculate on model
    ##   [Note] Any error reports that follow may simply reflect missing values in model variables.

    ## Checking model sizes and dimensions

    ##   [Note] This model is not fully initialized. This is not an error.
    ##          To see which variables are not initialized, use model$initializeInfo().
    ##          For more information on model initialization, see help(modelInitialization).

    ## ===== Monitors =====
    ## thin = 1: alpha, alpha_prior, Gamma, Gamma_prior, mu, mu_prior, omega, omega_prior, sigma, sigma_prior
    ## ===== Samplers =====
    ## RW sampler (4)
    ##   - sigma
    ##   - Gamma
    ##   - omega
    ##   - mu
    ## posterior_predictive_branch sampler (2)
    ##   - omega_prior
    ##   - alpha_prior
    ## posterior_predictive sampler (3)
    ##   - sigma_prior
    ##   - Gamma_prior
    ##   - mu_prior
    ## conjugate sampler (41)
    ##   - alpha
    ##   - beta[]  (10 elements)
    ##   - phi[]  (30 elements)
    ## binary sampler (10)
    ##   - delta[]  (10 elements)

    ## Defining model

    ## Building model

    ## Setting data and initial values

    ## Running calculate on model
    ##   [Note] Any error reports that follow may simply reflect missing values in model variables.

    ## Checking model sizes and dimensions

    ##   [Note] This model is not fully initialized. This is not an error.
    ##          To see which variables are not initialized, use model$initializeInfo().
    ##          For more information on model initialization, see help(modelInitialization).

    ## Checking model calculations

    ## NAs were detected in model variables: sigma_prior, logProb_sigma_prior, Gamma_prior, logProb_Gamma_prior, omega_prior, logProb_omega_prior, alpha_prior, logProb_alpha_prior, mu_prior, logProb_mu_prior, delta_prior, logProb_delta_prior, lifted_omega_prior_times__oP_oP1_minus_delta_prior_cP_times_0_dot_0001_plus_delta_prior_times_1000_cP, beta_prior, logProb_beta_prior.

    ## Compiling
    ##   [Note] This may take a minute.
    ##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.

    ## running chain 1...

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|

    ## running chain 2...

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|

    ## running chain 3...

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|
    ## There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.

``` r
posterior_on_number_of_covariates_selected(mcmc_out, warmup_prop = 0)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
traceplot(mcmc_out)
```

![](README_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
posterior_predictive_estimate(data = data, mcmc_out = mcmc_out)
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

### Prior posterior plot

``` r
PriorPosteriorPlotScalarParams(mcmc_out, parnames = c("sigma", "Gamma", "omega", "alpha", "mu")) + scale_x_log10()
```

![](README_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
PriorPosteriorPlotLargeVectorParams(mcmc_out, parnames = c("beta", "delta")) 
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
PriorPosteriorPlotLargeVectorParams(mcmc_out, parnames = c("beta")) + scale_x_log10()
```

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning in self$trans$transform(x): NaNs produced

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 3791 rows containing non-finite values (stat_density).

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 2 rows containing non-finite values (stat_bin).

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

## Example 2: regression on the delay parameter, after reparametrisation

``` r
set.seed(0)
data = SimulateDataOneVariableParam_reparam(beta = c(0.3, 0.3, 0.15, 0, 0, 0, 0, 0, 0), psi2 = 0.25, psi3 = 0.15, sigma = 1.)
```

    ## Joining, by = "i"

``` r
data
```

    ## # A tibble: 300 × 28
    ##    times     i phi_i  gphi     Y    V2     V3    V4    V5     V6    V7      V8
    ##    <dbl> <int> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>
    ##  1   0       1  1.39  1.34  1.22  1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  2  22.2     1  1.39  2.56  4.03  1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  3  44.4     1  1.39  4.52  5.20  1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  4  66.7     1  1.39  7.13  9.08  1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  5  88.9     1  1.39  9.83  9.56  1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  6 111.      1  1.39 12.0  10.7   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  7 133.      1  1.39 13.4  13.0   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  8 156.      1  1.39 14.2  14.3   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ##  9 178.      1  1.39 14.6  14.4   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ## 10 200       1  1.39 14.8  14.4   1.26 -0.236 0.359  1.30 -0.119  1.47 -0.0704
    ## # … with 290 more rows, and 16 more variables: V9 <dbl>, V10 <dbl>,
    ## #   beta1 <dbl>, beta2 <dbl>, beta3 <dbl>, beta4 <dbl>, beta5 <dbl>,
    ## #   beta6 <dbl>, beta7 <dbl>, beta8 <dbl>, beta9 <dbl>, Lmax <dbl>, Tmax <dbl>,
    ## #   psi1 <dbl>, psi2 <dbl>, psi3 <dbl>

### View simulated data

``` r
PlotSimulatedData(data)
```

![](README_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
mcmc_out = FitOneVariableParamSpikeAndSlabDelayReparam(data, nu_0 = 10^(-4), niter = 15000)
```

    ## Defining model

    ## Building model

    ## Setting data and initial values

    ## Running calculate on model
    ##   [Note] Any error reports that follow may simply reflect missing values in model variables.

    ## Checking model sizes and dimensions

    ## ===== Monitors =====
    ## thin = 1: alpha, alpha_prior, Gamma2, Gamma2_prior, psi1, psi1_prior, psi2, psi2_prior, psi3, psi3_prior, sigma2, sigma2_prior
    ## ===== Samplers =====
    ## RW sampler (35)
    ##   - psi1
    ##   - psi2
    ##   - psi3
    ##   - Gamma2
    ##   - sigma2
    ##   - phi[]  (30 elements)
    ## posterior_predictive sampler (6)
    ##   - sigma2_prior
    ##   - Gamma2_prior
    ##   - alpha_prior
    ##   - psi1_prior
    ##   - psi2_prior
    ##   - psi3_prior
    ## conjugate sampler (10)
    ##   - alpha
    ##   - beta[]  (9 elements)
    ## binary sampler (9)
    ##   - delta[]  (9 elements)

    ## Defining model

    ## Building model

    ## Setting data and initial values

    ## Running calculate on model
    ##   [Note] Any error reports that follow may simply reflect missing values in model variables.

    ## Checking model sizes and dimensions

    ## Checking model calculations

    ## Compiling
    ##   [Note] This may take a minute.
    ##   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.

    ## running chain 1...

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|

    ## running chain 2...

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|

    ## running chain 3...

    ## |-------------|-------------|-------------|-------------|
    ## |-------------------------------------------------------|
    ## There are individual pWAIC values that are greater than 0.4. This may indicate that the WAIC estimate is unstable (Vehtari et al., 2017), at least in cases without grouping of data nodes or multivariate data nodes.

``` r
posterior_on_number_of_covariates_selected(mcmc_out)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
traceplot(mcmc_out)
```

![](README_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
posterior_predictive_estimate(data = data, mcmc_out = mcmc_out)
```

![](README_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
PriorPosteriorPlotScalarParams(mcmc_out = mcmc_out, parnames = c("alpha", "sigma2", "psi1", "psi2", "psi3", "Gamma2"))
```

![](README_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
PriorPosteriorPlotScalarParams(mcmc_out = mcmc_out, parnames = c("alpha", "sigma2", "psi1", "psi2", "psi3", "Gamma2")) + scale_x_log10()
```

![](README_files/figure-markdown_github/unnamed-chunk-22-1.png)
