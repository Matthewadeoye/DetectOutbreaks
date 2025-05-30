% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{infer.old}
\alias{infer.old}
\title{A function to carryout Bayesian inference for the spatio-temporal models presented in (Knorr-Held, et.al., 2003).}
\usage{
infer.old(
  y,
  e_it,
  Model,
  adjmat,
  num_iteration = 30000,
  Stan = TRUE,
  GPU = FALSE,
  nchains = 4,
  iter = 4000,
  seed = NULL,
  verbose = F,
  ModEvid = F,
  OutbreakProb = F,
  adaptdelta = 0.9,
  Burn.in = 1000
)
}
\arguments{
\item{y}{A space-time data matrix (locations on the rows, and time on the columns).}

\item{e_it}{A space-time data matrix showing the number of susceptible individuals (locations on the rows, and time on the columns).}

\item{Model}{The model specification (ranges from 0 to 7).}

\item{adjmat}{The adjacency matrix describing the connectivity/neighborhood structure of the spatial locations.}

\item{num_iteration}{The required number of iterations for inference via MCMC. Default set to 30000.}

\item{Stan}{A logical argument asking whether to fit the model using Stan's HMC via cmdstanr.}

\item{GPU}{A logical argument for accelerating Stan's computations using OpenCL (works with GPU only).}

\item{nchains}{The number of parallel MCMC chains to run using Stan's HMC.}

\item{iter}{The required number of iterations for inference via HMC.}

\item{seed}{A random number for reproducibility.}

\item{verbose}{A logical argument asking whether to display (unnecessary) outputs.}

\item{ModEvid}{A logical argument asking whether to compute and print the model evidence.}

\item{OutbreakProb}{A logical argument asking whether to compute and return the posterior marginal probabilities of outbreak.}

\item{adaptdelta}{Target average proposal acceptance probability during Stan's adaptation period. Default set to 0.90.}

\item{Burn.in}{MCMC burn-in period. Default set to 3000.}
}
\value{
Samples from the target posterior distribution (and outbreak probabilities, model evidence).
}
\description{
A function to carryout Bayesian inference for the spatio-temporal models presented in (Knorr-Held, et.al., 2003).
}
\examples{
set.seed(4);
sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
TRUTHS1<- simulate.old(Model = 0, time = 48, adj.matrix = sim_adjmat);
HMCfit1<-infer.old(TRUTHS1[[1]],TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=TRUE);
HMCfit1<- HMCfit1[[1]]
mcmc.plot(HMCfit1)
inf.plot(HMCfit1)

MCMCfit1<-infer.old(TRUTHS1[[1]],TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,Stan=FALSE,ModEvid=TRUE,OutbreakProb=TRUE);

}
