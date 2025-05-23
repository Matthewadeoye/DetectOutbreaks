% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OutbreakProbability.R
\name{Decoding}
\alias{Decoding}
\title{A function to compute the posterior marginal probabilities of outbreak.}
\usage{
Decoding(y, e_it, r, s, u, Gamma, B, Model, adjmat, Cyclic = F)
}
\arguments{
\item{y}{A space-time data matrix (locations on the row, and time on the columns).}

\item{e_it}{A space-time data matrix showing the number of susceptible individuals.}

\item{r}{Trend component.}

\item{s}{Seasonal component.}

\item{u}{Spatial component.}

\item{Gamma}{Transition probability matrix.}

\item{B}{Autoregressive coefficients.}

\item{Model}{The model specification (ranges from 0 to 6).}

\item{adjmat}{The adjacency matrix describing the connectivity of spatial locations.}

\item{Cyclic}{A logical argument asking whether a cyclic RW1 prior on s_t was used for inference.}
}
\value{
A matrix of outbreak probabilities.
}
\description{
A function to compute the posterior marginal probabilities of outbreak.
}
\examples{
set.seed(4);
sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
TRUTHS<- simulate(Model = 1, time = 48, adj.matrix = sim_adjmat, T.prob = G(0.2, 0.4), B = 0.65)
singleOutProb<- Decoding(y=TRUTHS[[1]],e_it=TRUTHS[[2]],r=TRUTHS[[3]],s=TRUTHS[[4]],u=TRUTHS[[5]],Gamma=G(0.2, 0.4),B=0.65,Model=1,Cyclic=TRUE)
image(t(singleOutProb))

set.seed(4);
TRUTHS1<- simulate.old(Model = 1, time = 48, adj.matrix = sim_adjmat, T.prob = G(0.2, 0.4), B = 0.65)
singleOutProb<-Decoding(y=TRUTHS1[[1]],e_it=TRUTHS1[[2]],r=TRUTHS1[[3]],s=TRUTHS1[[4]],u=TRUTHS1[[5]],Gamma=G(0.2,0.4),B=0.65,Model=1,Cyclic=FALSE)
image(t(singleOutProb))

}
