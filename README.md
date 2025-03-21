
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DetectOutbreaks

<!-- badges: start -->
<!-- badges: end -->

The goal of DetectOutbreaks is to introduce a novel approach to the
analysis of a spatio-temporal model designed for disease outbreak
surveillance. The package is built to help facilitate model simulations
and efficient Bayesian inference methods.

# Dependencies

<!-- badges: start -->
<!-- badges: end -->

DetectOutbreaks depends on Stan through its R interface, cmdstanr. If
not already installed, you can do so using

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
install_cmdstan()
```

## Installation

You can install the development version of DetectOutbreaks from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Matthewadeoye/DetectOutbreaks")
```

## Simulation, visualization, inference, and outbreak detection.

This is a basic example which shows you how to simulate and visualize
space-time data from the null model:

``` r
library(DetectOutbreaks)
#> 
#> Attaching package: 'DetectOutbreaks'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
set.seed(0); TRUTHS<- simulate(Model = 0, time = 60, adj.matrix = sim_adjmat)
sim.plot(sim.object = TRUTHS)
```

<img src="man/figures/README-simulation-1.png" width="100%" /><img src="man/figures/README-simulation-2.png" width="100%" />

<img src="man/figures/README-figures-1.png" width="100%" /><img src="man/figures/README-figures-2.png" width="100%" />

We can perform Bayesian inference using either HMC in Stan or a bespoke
MCMC algorithm in DetectOutbreaks:

``` r
#HMCfit<- infer(y = TRUTHS[[1]], e_it = TRUTHS[[2]], Model = 0, adjmat = sim_adjmat, Stan = TRUE, ModEvid = TRUE, OutbreakProb = FALSE)
#mcmc.plot(HMCfit)
#inf.plot(HMCfit)

#MCMCfit<- infer(y = TRUTHS[[1]], e_it = TRUTHS[[2]], Model = 0, adjmat = sim_adjmat, Stan = FALSE, ModEvid = TRUE, OutbreakProb = FALSE)
#mcmc.plot(MCMCfit)
#inf.plot(MCMCfit)
```

We can also perform simulation, inference and detect outbreaks in the
simulated datasets using more complex models in DetectOutbreaks:

``` r
library(DetectOutbreaks)
sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
set.seed(0); TRUTHS1<- simulate(Model = 1, time = 60, adj.matrix = sim_adjmat)
#sim.plot(sim.object = TRUTHS1)
#HMCfit1<- infer(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], Model = 1, adjmat = sim_adjmat, Stan = TRUE, ModEvid = TRUE, OutbreakProb = FALSE)
#mcmc.plot(HMCfit1)
#inf.plot(HMCfit1)
#MarginalProbabilities<- OutbreakProbability(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], inf.object = HMCfit1, adjmat = sim_adjmat, Model = 1)
#image(t(MarginalProbabilities))

#MCMCfit1<- infer(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], Model = 1, adjmat = sim_adjmat, Stan = FALSE, ModEvid = TRUE, OutbreakProb = FALSE)
#mcmc.plot(MCMCfit1)
#inf.plot(MCMCfit1)
#MarginalProbabilities<- OutbreakProbability(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], inf.object = #MCMCfit1, adjmat = sim_adjmat, Model = 1)
#image(t(MarginalProbabilities))
```

We can compute model evidence (also known as marginal log-likelihood)
through importance sampling Monte Carlo simulations for any model fitted
by DetectOutbreaks:

``` r
library(DetectOutbreaks)
sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
set.seed(0); TRUTHS1<- simulate(Model = 1, time = 60, adj.matrix = sim_adjmat)
#HMCfit1<- infer(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], Model = 1, adjmat = sim_adjmat, Stan = TRUE, ModEvid = FALSE, OutbreakProb = FALSE)
#ModelEvidence(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], adjmat = sim_adjmat, Model = 1, inf.object = HMCfit1, num_samples = 50000)

#MCMCfit1<- infer(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], Model = 1, adjmat = sim_adjmat, Stan = FALSE, ModEvid = FALSE, OutbreakProb = FALSE)
#ModelEvidence(y = TRUTHS1[[1]], e_it = TRUTHS1[[2]], adjmat = sim_adjmat, Model = 1, inf.object = MCMCfit1, num_samples = 50000)
```
