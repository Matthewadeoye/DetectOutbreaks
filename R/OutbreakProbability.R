
#FORWARD FILTER for old model
forwardfilter <- function(y, e_it, r, s, u, Gamma, B, Model, adjmat) {

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  ndept <- nrow(y)
  nstate <- ncol(Gamma)
  time <- ncol(y)
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init_density<- state_dist_cpp(Gamma[1, 2], Gamma[2, 1])
  init_density<- log(init_density)

  y<- ifelse(is.na(y), -1, y)

  if(Model %in% c(0,1,2,4,5,7)){
    AllForwardprobs<- vector("list", ndept)

    for (i in 1:ndept) {
      Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
      if(y[i, 1] == -1){
        alpha.1 <- init_density[1]
        alpha.2 <- init_density[2]
      }else{
        alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i, 1]), log = TRUE)
      }
      Forwardprob[1, ] <- c(alpha.1, alpha.2)

      for (t in 2:time) {
        if(y[i, t] == -1){
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }else{
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i,t]), log = TRUE)))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }
      }
      AllForwardprobs[[i]]<- Forwardprob
    }

    return(AllForwardprobs)
  }

  else if(Model %in% c(3,6)){
    AllForwardprobs<- vector("list", ndept)

    for (i in 1:ndept) {
      Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
      if(y[i, 1] == -1){
        alpha.1 <- init_density[1]
        alpha.2 <- init_density[2]
      }else{
      alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
      alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i,1] + B[2] * z_it2[i,1]), log = TRUE)
    }
      Forwardprob[1, ] <- c(alpha.1, alpha.2)

      for (t in 2:time) {
        if(y[i, t] == -1){
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }else{
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i,t] + B[2] * z_it2[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[t] + u[i] + B[1] * z_it[i,t] + B[2] * z_it2[i,t]), log = TRUE)))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }
      }
      AllForwardprobs[[i]]<- Forwardprob
    }

    return(AllForwardprobs)
  }
}

#BACKWARD SWEEP for old model
backwardsweep <- function(y, e_it, r, s, u, Gamma, B, Model, adjmat) {

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  ndept <- nrow(y)
  nstate <- ncol(Gamma)
  time <- ncol(y)
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])

  y<- ifelse(is.na(y), -1, y)

  if(Model %in% c(0,1,2,4,5,7)){
    Allbackwardprob<- vector("list", ndept)

    for (i in 1:ndept) {
      Backwardprob<- matrix(NA, nrow = time, ncol = nstate)
      beta.1 <- 0
      beta.2 <- 0
      Backwardprob[time, ] <- c(beta.1, beta.2)

      for (t in (time-1):1) {
        if(y[i, t+1] == -1){
          Betas<- c(beta.1, beta.2)
          beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11))
          beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21))
          Backwardprob[t, ] <- c(beta.1, beta.2)
        }else{
          Betas<- c(beta.1, beta.2)
          beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]), log = TRUE), Betas[2] + gamma_12 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]), log = TRUE)))
          beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + B[1] * z_it[i, t+1]), log = TRUE), Betas[2] + gamma_22 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + B[1] * z_it[i, t+1]), log = TRUE)))
          Backwardprob[t, ] <- c(beta.1, beta.2)
        }
      }
      Allbackwardprob[[i]] <- Backwardprob
    }
    return(Allbackwardprob)
  }

  else if(Model %in% c(3,6)){
    Allbackwardprob<- vector("list", ndept)

    for (i in 1:ndept) {
      Backwardprob<- matrix(NA, nrow = time, ncol = nstate)
      beta.1 <- 0
      beta.2 <- 0
      Backwardprob[time, ] <- c(beta.1, beta.2)

      for (t in (time-1):1) {
        if(y[i, t+1] == -1){
          Betas<- c(beta.1, beta.2)
          beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11))
          beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21))
          Backwardprob[t, ] <- c(beta.1, beta.2)
        }else{
          Betas<- c(beta.1, beta.2)
          beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]), log = TRUE), Betas[2] + gamma_12 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i]), log = TRUE)))
          beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + B[1] * z_it[i,t+1] + B[2] * z_it2[i,t+1]), log = TRUE), Betas[2] + gamma_22 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[t+1] + u[i] + B[1] * z_it[i,t+1] + B[2] * z_it2[i,t+1]), log = TRUE)))
          Backwardprob[t, ] <- c(beta.1, beta.2)
        }
      }
      Allbackwardprob[[i]] <- Backwardprob
    }
    return(Allbackwardprob)
  }
}

#FORWARD FILTER for new (cyclic) model
forwardfilter2<- function(y, e_it, r, s, u, Gamma, B, Model, adjmat) {

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  ndept <- nrow(y)
  nstate <- ncol(Gamma)
  time <- ncol(y)
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init_density<- state_dist_cpp(Gamma[1, 2], Gamma[2, 1])
  init_density<- log(init_density)

  y<- ifelse(is.na(y), -1, y)

  if(Model %in% c(0,1,2,4,5,7)){
    AllForwardprobs<- vector("list", ndept)

    for (i in 1:ndept) {
      Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
      if(y[i, 1] == -1){
        alpha.1 <- init_density[1]
        alpha.2 <- init_density[2]
      }else{
      alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
      alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i, 1]), log = TRUE)
}
      Forwardprob[1, ] <- c(alpha.1, alpha.2)

      for (t in 2:time) {
        if(y[i, t] == -1){
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }else{
        month_index<- (t - 1) %% 12 + 1
        Alphas<- c(alpha.1, alpha.2)
        alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE)))
        alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t]), log = TRUE)))
        Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }
      }
      AllForwardprobs[[i]]<- Forwardprob
    }

    return(AllForwardprobs)
  }

  else if(Model %in% c(3,6)){
    AllForwardprobs<- vector("list", ndept)

    for (i in 1:ndept) {
      Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
      if(y[i, 1] == -1){
        alpha.1 <- init_density[1]
        alpha.2 <- init_density[2]
      }else{
      alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
      alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i,1] + B[2] * z_it2[i,1]), log = TRUE)
      }

      Forwardprob[1, ] <- c(alpha.1, alpha.2)

      for (t in 2:time) {
        if(y[i, t] == -1){
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }else{
        month_index<- (t - 1) %% 12 + 1
        Alphas<- c(alpha.1, alpha.2)
        alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE)))
        alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t] + B[2] * z_it2[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t] + B[2] * z_it2[i,t]), log = TRUE)))
        Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }
      }
      AllForwardprobs[[i]]<- Forwardprob
    }

    return(AllForwardprobs)
  }
}

#BACKWARD SWEEP for new (cyclic) model
backwardsweep2<- function(y, e_it, r, s, u, Gamma, B, Model, adjmat) {

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  ndept <- nrow(y)
  nstate <- ncol(Gamma)
  time <- ncol(y)
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])

  y<- ifelse(is.na(y), -1, y)

  if(Model %in% c(0,1,2,4,5,7)){
    Allbackwardprob<- vector("list", ndept)

    for (i in 1:ndept) {
      Backwardprob<- matrix(NA, nrow = time, ncol = nstate)
      beta.1 <- 0
      beta.2 <- 0
      Backwardprob[time, ] <- c(beta.1, beta.2)

      for (t in (time-1):1) {
        if(y[i, t+1] == -1){
          Betas<- c(beta.1, beta.2)
          beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11))
          beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21))
          Backwardprob[t, ] <- c(beta.1, beta.2)
        }else{
        month_index<- t %% 12 + 1
        Betas<- c(beta.1, beta.2)
        beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i]), log = TRUE), Betas[2] + gamma_12 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i]), log = TRUE)))
        beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i] + B[1] * z_it[i, t+1]), log = TRUE), Betas[2] + gamma_22 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i] + B[1] * z_it[i, t+1]), log = TRUE)))
        Backwardprob[t, ] <- c(beta.1, beta.2)
        }
      }
      Allbackwardprob[[i]] <- Backwardprob
    }
    return(Allbackwardprob)
  }

  else if(Model %in% c(3,6)){
    Allbackwardprob<- vector("list", ndept)

    for (i in 1:ndept) {
      Backwardprob<- matrix(NA, nrow = time, ncol = nstate)
      beta.1 <- 0
      beta.2 <- 0
      Backwardprob[time, ] <- c(beta.1, beta.2)

      for (t in (time-1):1) {
        if(y[i, t+1] == -1){
          Betas<- c(beta.1, beta.2)
          beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11))
          beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21))
          Backwardprob[t, ] <- c(beta.1, beta.2)
        }else{
        month_index<- t %% 12 + 1
        Betas<- c(beta.1, beta.2)
        beta.1 <- logSumExp_cpp(c(Betas[1] + gamma_11 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i]), log = TRUE), Betas[2] + gamma_12 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i]), log = TRUE)))
        beta.2 <- logSumExp_cpp(c(Betas[1] + gamma_21 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i] + B[1] * z_it[i,t+1] + B[2] * z_it2[i,t+1]), log = TRUE), Betas[2] + gamma_22 + dpois(y[i, t+1], lambda = e_it[i, t+1] * exp(r[t+1] + s[month_index] + u[i] + B[1] * z_it[i,t+1] + B[2] * z_it2[i,t+1]), log = TRUE)))
        Backwardprob[t, ] <- c(beta.1, beta.2)
        }
      }
      Allbackwardprob[[i]] <- Backwardprob
    }
    return(Allbackwardprob)
  }
}

#' A function to compute the posterior marginal probabilities of outbreak.
#'
#' @param y A space-time data matrix (locations on the row, and time on the columns).
#' @param e_it A space-time data matrix showing the number of susceptible individuals.
#' @param r Trend component.
#' @param s Seasonal component.
#' @param u Spatial component.
#' @param Gamma Transition probability matrix.
#' @param B Autoregressive coefficients.
#' @param Model The model specification (ranges from 0 to 6).
#' @param Cyclic A logical argument asking whether a cyclic RW1 prior on s_t was used for inference.
#' @param adjmat The adjacency matrix describing the connectivity of spatial locations.
#'
#' @return A matrix of outbreak probabilities.
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 1, time = 48, adj.matrix = sim_adjmat, T.prob = G(0.2, 0.4), B = 0.65)
#' singleOutProb<- Decoding(y=TRUTHS[[1]],e_it=TRUTHS[[2]],r=TRUTHS[[3]],s=TRUTHS[[4]],u=TRUTHS[[5]],Gamma=G(0.2, 0.4),B=0.65,Model=1,Cyclic=TRUE)
#' image(t(singleOutProb))
#'
#' set.seed(4);
#' TRUTHS1<- simulate.old(Model = 1, time = 48, adj.matrix = sim_adjmat, T.prob = G(0.2, 0.4), B = 0.65)
#' singleOutProb<-Decoding(y=TRUTHS1[[1]],e_it=TRUTHS1[[2]],r=TRUTHS1[[3]],s=TRUTHS1[[4]],u=TRUTHS1[[5]],Gamma=G(0.2,0.4),B=0.65,Model=1,Cyclic=FALSE)
#' image(t(singleOutProb))
#'
Decoding <- function(y, e_it, r, s, u, Gamma, B, Model, adjmat, Cyclic = F) {
  ndept<- nrow(y)
  time <- ncol(y)
  if(!Cyclic){
  Allforwardprobs<- forwardfilter(y, e_it, r, s, u, Gamma, B, Model, adjmat)
  Allbackwardprobs<- backwardsweep(y, e_it, r, s, u, Gamma, B, Model, adjmat)
  Res<- matrix(NA, ndept, time)
  for(i in 1:ndept){
    for(j in 1:time){
      P1<- Allforwardprobs[[i]][j,1] + Allbackwardprobs[[i]][j,1]
      P2<- Allforwardprobs[[i]][j,2] + Allbackwardprobs[[i]][j,2]
      Res[i,j]<- exp(P2 - logSumExp_cpp(c(P1,P2)))
      #browser()
    }
  }
  return(Res)
  }else{
    Allforwardprobs<- forwardfilter2(y, e_it, r, s, u, Gamma, B, Model, adjmat)
    Allbackwardprobs<- backwardsweep2(y, e_it, r, s, u, Gamma, B, Model, adjmat)
    Res<- matrix(NA, ndept, time)
    for(i in 1:ndept){
      for(j in 1:time){
        P1<- Allforwardprobs[[i]][j,1] + Allbackwardprobs[[i]][j,1]
        P2<- Allforwardprobs[[i]][j,2] + Allbackwardprobs[[i]][j,2]
        Res[i,j]<- exp(P2 - logSumExp_cpp(c(P1,P2)))
        #browser()
      }
    }
    return(Res)
 }
}


#' Compute posterior marginal probabilities of outbreak using all thinned MCMC samples.
#'
#' @param y A space-time data matrix (locations on the row, and time on the columns).
#' @param e_it A space-time data matrix showing the number of susceptible individuals.
#' @param inf.object Resulting object from either "infer" or "infer.old" functions.
#' @param adjmat The adjacency matrix describing the connectivity of spatial locations.
#' @param Model The model specification (ranges from 0 to 6).
#' @param burn.in Burn-in period for MCMC.
#' @param Cyclic A logical argument asking whether a cyclic RW1 prior on s_t was used for inference.
#'
#' @return A matrix of outbreak probabilities.
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 1, time = 48, adj.matrix = sim_adjmat, T.prob = G(0.2, 0.4), B = 0.65)
#' fit<-infer(TRUTHS[[1]],TRUTHS[[2]],Model=1,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=FALSE)
#' OutP<-OutbreakProbability(y=TRUTHS[[1]],e_it=TRUTHS[[2]],inf.object=fit,adjmat=sim_adjmat,Model=1,Cyclic=TRUE)
#' image(t(OutP))
#'
OutbreakProbability<- function(y, e_it, inf.object, adjmat, Model, burn.in = 2000, Cyclic = T){
  time<- ncol(y)
  ndept<- nrow(y)

  y<- ifelse(is.na(y), -1, y)

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  if(!is.data.frame(inf.object)){
    fullG12.draws<- stack(as.data.frame(inf.object$draws(variables = "G12")[,1,]))[,1]
    fullG21.draws<- stack(as.data.frame(inf.object$draws(variables = "G21")[,1,]))[,1]
    fullr.draws<- as.data.frame(inf.object$draws(variables = "r")[,1,])
    fulls.draws<- as.data.frame(inf.object$draws(variables = "s")[,1,])
    fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
    fulld.draws<- stack(as.data.frame(inf.object$draws(variables = "state1_stationary_dist")[,1,]))[,1]
  }else{
    fullG12.draws<- as.numeric(inf.object[-(1:burn.in), 1])
    fullG21.draws<- as.numeric(inf.object[-(1:burn.in), 2])
    fullr.draws<- inf.object[-(1:burn.in), 5+(1:time)]
    fulls.draws<- inf.object[-(1:burn.in), 5+time+(1:12)]
    fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
    fulld.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+2+1]
  }

  thinning<- numeric(floor(nrow(fullr.draws)/10))
  thinning[1]<- 10
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + 10
  }

  G12.draws<- fullG12.draws[thinning]
  G21.draws<- fullG21.draws[thinning]
  r.draws<- fullr.draws[thinning, ]
  s.draws<- fulls.draws[thinning, ]
  u.draws<- fullu.draws[thinning, ]
  d.draws<- fulld.draws[thinning]

  Ex_Xit<- matrix(0, nrow = ndept, ncol = time)

 if(Model %in% c(1,2,4,5,7)){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    for(index in 1:length(thinning)){
      r<- as.numeric(r.draws[index,])
      s<- as.numeric(s.draws[index,])
      u<- as.numeric(u.draws[index,])
      Decode_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B.draws[index], Model = Model, adjmat = adjmat, Cyclic = Cyclic)
      Ex_Xit<- Ex_Xit + Decode_Xit
    }
  }else if(Model %in% c(3,6)){
    if(!is.data.frame(inf.object)){
      B.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
    }else{
      B.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+(1:2)]
    }
    B.draws<- B.draws[thinning, ]
    for(index in 1:length(thinning)){
      r<- as.numeric(r.draws[index,])
      s<- as.numeric(s.draws[index,])
      u<- as.numeric(u.draws[index,])
      B<- as.numeric(B.draws[index,])
      Decode_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B, Model = Model, adjmat = adjmat, Cyclic = Cyclic)
      Ex_Xit<- Ex_Xit + Decode_Xit
    }
  }
  Ex_Xit<- Ex_Xit/length(thinning)
  return(Ex_Xit)
}
