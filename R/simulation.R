
sim.AR2 <- function(time, ar1, ar2, sd, init.r1 = -14, init.r2 = -14){
  rt<- numeric(time)
  rt[1] <- init.r1
  rt[2] <- init.r2
  for (i in 2:(time-1)){
    rt[i+1] <- ar1*rt[i] + ar2*rt[i-1] + rnorm(1, mean=0, sd=sd)
  }
  return(rt)
}

sim.RW2<- function(time, sd=0.01, init.r1 = -14, init.r2 = -14){
  r <- numeric(time)

  r[1] <- init.r1
  r[2] <- init.r2

  for(t in 3:time){
    epsilon <- rnorm(1, mean = 0, sd = sd)
    r[t] <- 2*(r[t - 1]) - r[t - 2] + epsilon
  }
  return(r)
}

sim.Seasonals<- function(time,  sd = 0.03){
  s <- numeric(time)
  initial_WN <- rnorm(time, mean = 0, sd = sd/time)

  for (i in 1:time) {
    if (i <= 11) {
      s[i] <- initial_WN[i]
    }else {
      s[i] <- rnorm(1, mean = -sum(s[(i-1):(i-11)]), sd = sd/time)
    }
  }
  return(s)
}

sim.Seasonals2<- function(Amplitude, Cycle = 1:12){
  frequency <- 1/length(Cycle)
  s <- sin(2 * pi * frequency * Cycle)
  s <- Amplitude * s
  return(s)
}

sim.GMRF <- function(n, Q, tol=1e-12){
  Es <- eigen(Q)
  LAMBDAS <- rev(Es$values)
  ai <- LAMBDAS > tol
  Vtilde <- Es$vectors[,sum(ai):1]

  yij <- matrix(
    stats::rnorm(sum(ai) * n, 0, sqrt(rep(LAMBDAS[ai], n)^-1)),
    nrow=sum(ai),
    ncol=n)

  t(Vtilde %*% yij)
}

sim.Seasonals3<- function(sd = 0.08){

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), 2,-1)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW1PrecMat<- (1/sd^2) * RW1PrecMat
  s<- c(sim.GMRF(1, RW1PrecMat))
  return(s)
}

sim.Spatials<- function(adj.matrix, sd = 0.2){
  Q<- -1 * adj.matrix
  diag(Q)<- -rowSums(Q, na.rm = T)
  Q<- (1/sd^2)*Q
  u<- c(sim.GMRF(1,Q))
  return(u)
}

MChain<- function(time, G12, G21){

  transition_matrix <- matrix(c(1-G12, G12, G21, 1-G21), nrow = 2, byrow = TRUE)

  initial_state <- 0

  states <- numeric(time)

  states[1] <- initial_state
  for(i in 2:time) {
    current_state <- states[i - 1]
    next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
    states[i] <- next_state
  }
  return(states)
}


#' A function to simulate spatio-temporal datasets from the spatio-temporal models presented in (Adeoye, et.al., 2024).
#'
#' @param Model The model specification (ranges from 0 to 7).
#' @param time Represent the time points.
#' @param adj.matrix The adjacency matrix describing the connectivity of spatial locations.
#' @param e_it A space-time data matrix showing the number of susceptible individuals.
#' @param B Autoregressive coefficient(s).
#' @param T.prob Transition probability matrix.
#' @param r Trend component.
#' @param s Seasonal component.
#' @param u Spatial component.
#'
#' @return Spatio-temporal dataset
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 0, time = 48, adj.matrix = sim_adjmat)
#'
simulate<- function(Model, time, adj.matrix,
                     e_it=matrix(c(rep(c(rpois(time, 500000), rpois(time, 1000000)), 4), rpois(time, 500000)),
                                 byrow = T, ncol = time),
                     B = c(1.68, 0.20), T.prob = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T),
                     r = sim.RW2(time, sd=0.009), s = sim.Seasonals2(Amplitude = 1.4),
                     u = sim.Spatials(adj.matrix)){
  #sim.AR2(time, ar1=0.95, ar2=0.065, sd=0.57)
  #Models 1,2,7 B = 1.65
  #Models 2 B = 1.68 for  model-evidence study
  #Models 4 B = 0.55
  #Models 5 B = 0.45
  #Model 3 B = c(1.25, 0.75)
  #Model 6 B=  c(0.35, 0.20)
  ndept<- nrow(adj.matrix)
  y_it<- matrix(NA, ndept, time)
  EpidemicIndicator<- matrix(NA, ndept, time)

  if(Model == 0){
    for(i in 1:ndept){
      for(t in 1:time){
        m<- (t - 1) %% 12 + 1
        lograte <- r[t] + s[m] + u[i]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 1){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        z_it<- ifelse(y_it[i, t-1]>0, 1, 0)
        lograte<- r[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 2){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it<-  ifelse(y_it[i, t-1]>0 || any(y_it[indexes, t-1]>0), 1, 0)
        lograte<- r[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 3){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it1<- ifelse(y_it[i, t-1]>0, 1, 0)
        z_it2<- ifelse(any(y_it[indexes, t-1] > 0), 1, 0)
        lograte <- r[t] + s[m] + u[i] + (z_it1 * EpidemicIndicator[i, t] * B[1]) + (z_it2 * EpidemicIndicator[i, t] * B[2])
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 4){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        z_it<- log(y_it[i, t-1] + 1)
        lograte<- r[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 5){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it<- log(y_it[i, t-1] + sum(y_it[indexes, t-1]) + 1)
        lograte<- r[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 6){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it1<- log(y_it[i, t-1] + 1)
        z_it2<- log(sum(y_it[indexes, t-1]) + 1)
        lograte <- r[t] + s[m] + u[i] + (z_it1 * EpidemicIndicator[i, t] * B[1]) + (z_it2 * EpidemicIndicator[i, t] * B[2])
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 7){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        z_it<- 1
        lograte<- r[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }
}



#' A function to simulate spatio-temporal datasets from the spatio-temporal models presented in (Knorr-Held, et.al., 2003).
#'
#' @param Model The model specification (ranges from 0 to 7).
#' @param time Represent the time points
#' @param adj.matrix The adjacency matrix describing the connectivity of spatial locations.
#' @param e_it A space-time data matrix showing the number of susceptible individuals.
#' @param B Autoregressive coefficients.
#' @param T.prob Transition probability matrix.
#' @param r Trend component.
#' @param s Seasonal component.
#' @param u Spatial component.
#'
#' @return The simulated dataset.
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS1<- simulate.old(Model = 0, time = 48, adj.matrix = sim_adjmat)
#'
simulate.old<- function(Model, time, adj.matrix,
                    e_it=matrix(c(rep(c(rpois(time, 500000), rpois(time, 1000000)), 4), rpois(time, 500000)),
                                  byrow = T, ncol = time),
                    B = c(0.55, 0.45), T.prob = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T),
                    r = sim.AR2(time, 0.95, 0.065, 0.57), s = sim.Seasonals(time),
                    u = sim.Spatials(adj.matrix)){

  ndept<- nrow(adj.matrix)

  y_it<- matrix(NA, ndept, time)
  EpidemicIndicator<- matrix(NA, ndept, time)

  if(Model == 0){
    for(i in 1:ndept){
      for(t in 1:time){
        lograte <- r[t] + s[t] + u[i]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u))
  }

  else if(Model == 1){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        z_it<- ifelse(y_it[i, t-1]>0, 1, 0)
        lograte<- r[t] + s[t] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 2){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it<-  ifelse(y_it[i, t-1]>0 || any(y_it[indexes, t-1]>0), 1, 0)
        lograte<- r[t] + s[t] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 3){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it1<- ifelse(y_it[i, t-1]>0, 1, 0)
        z_it2<- ifelse(any(y_it[indexes, t-1] > 0), 1, 0)
        lograte <- r[t] + s[t] + u[i] + z_it1 * EpidemicIndicator[i, t] * B[1]
        + z_it2 * EpidemicIndicator[i, t] * B[2]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 4){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        z_it<- log(y_it[i, t-1] + 1)
        lograte<- r[t] + s[t] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 5){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it<- log(y_it[i, t-1] + sum(y_it[indexes, t-1]) + 1)
        lograte<- r[t] + s[t] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 6){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it1<- log(y_it[i, t-1] + 1)
        z_it2<- log(sum(y_it[indexes, t-1]) + 1)
        lograte <- r[t] + s[t] + u[i] + z_it1 * EpidemicIndicator[i, t] * B[1]
        + z_it2 * EpidemicIndicator[i, t] * B[2]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }

  else if(Model == 7){
    for(i in 1:ndept){
      EpidemicIndicator[i, ]<- MChain(time, T.prob[1, 2], T.prob[2, 1])
      y_it[i, 1]<- rpois(1, lambda = 1)
    }

    for(t in 2:time){
      for(i in 1:ndept){
        z_it<- 1
        lograte<- r[t] + s[t] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    return(list(y_it, e_it, r, s, u, EpidemicIndicator))
  }
}

