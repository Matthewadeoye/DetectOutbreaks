#Check valid simulations
check_sim <- function(sim.object, adjmat, Model) {
  y<- sim.object[[1]]
  e_it<- sim.object[[2]]
  r<- sim.object[[3]]
  s<- sim.object[[4]]
  u<- sim.object[[5]]
  x<- sim.object[[6]]
  ndept<- nrow(y)
  time<- ncol(y)

  lambda_it<- matrix(NA, nrow = ndept, ncol = time)
  S_naught<- matrix(NA, nrow = ndept, ncol = time)
  S_one<- matrix(NA, nrow = ndept, ncol = time)

  for(i in 1:ndept){
    for(t in 1:time){
      m<- (t - 1) %% 12 + 1
      lambda_it[i, t] = exp(r[t] + s[m] + u[i])
    }
  }

  if(Model != 0){
    for(i in 1:ndept){
      for(t in 1:time){
        if(x[i, t] == 0){
        S_naught[i, t]<- y[i, t] - e_it[i, t] * lambda_it[i, t]
        }else{
        S_one[i, t]<- y[i, t] - e_it[i, t] * lambda_it[i, t]
      }
    }
  }
  S_naught <- as.numeric(S_naught[!is.na(S_naught)])
  S_one<- as.numeric(S_one[!is.na(S_one)])

  vioSim<- data.frame(value = c(S_naught, S_one), group = factor(c(rep("S0", length(S_naught)), rep("S1", length(S_one)))))
  vioSim$group <- factor(vioSim$group, levels = unique(vioSim$group))
  Means<- c(mean(S_naught), mean(S_one))
  library(RColorBrewer)
  print(ggplot(vioSim, aes(x = group, y = value, fill = group)) +
          geom_violin(trim = FALSE, alpha = 0.7) +
          geom_point(x = 1, y = Means[1], color = "red", size = 3, shape = 18) +
          geom_point(x = 2, y = Means[2], color = "red", size = 3, shape = 18) +
          labs(title = "", x = "", y = "Value", fill = "Mean value") +
          theme_minimal() +
          scale_fill_brewer(palette = "Set2") +
          scale_x_discrete(labels = c("S0" = expression(S[O]),
                                      "S1" = expression(S[1]))) +
          scale_fill_manual(values = brewer.pal(5, "Set2"),
                            labels = c(expression(S[O]),
                                       expression(S[1])),
                            name = "Mean value"))  # Legend title
  }else{
    for(i in 1:ndept){
      for(t in 1:time){
          S_naught[i, t]<- y[i, t] - e_it[i, t] * lambda_it[i, t]
        }
    }
    S_naught <- as.numeric(S_naught[!is.na(S_naught)])
    vioSim<- data.frame(value = S_naught, group = factor(rep("S0", length(S_naught))))
    Means<- mean(S_naught)
    library(RColorBrewer)
    print(ggplot(vioSim, aes(x = group, y = value, fill = group)) +
            geom_violin(trim = FALSE, alpha = 0.7) +
            geom_point(x = 1, y = Means, color = "red", size = 3, shape = 18) +
            labs(title = "", x = "", y = "Value", fill = "Mean value") +
            theme_minimal() +
            scale_fill_brewer(palette = "Set2") +
            scale_x_discrete(labels = c("S0" = expression(S[O]))) +
            scale_fill_manual(values = brewer.pal(5, "Set2"),
                              labels = c(expression(S[O])),
                              name = "Mean value"))  # Legend title
  }
}


#Cholesky confirmation for using mvnfast
check_cholesky <- function(matrix) {
  result <- tryCatch({
    chol(matrix)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(result)
}

#Prior density for Trend components (r_t)
randomwalk2<- function(componentR, PrecisionR){
  time<- length(componentR)
  Sumres<- 0
  for(i in 3:time){
    res<- (componentR[i-2] - (2 * componentR[i-1]) + componentR[i])^2
    Sumres<- Sumres + res
  }
  return((time - 2)/2 * log(PrecisionR) - PrecisionR/2 * Sumres)
}

#Prior density for Seasonal components (s_t)
seasonalComp<- function(x, z){
  time<- length(x)
  Sumres<- 0
  for(i in 12:time){
    res<- (sum(x[(i-11):(i-0)]))^2
    Sumres<- Sumres + res
  }
  return((time - 11)/2 * log(z) - z/2 * Sumres)
}

#Cyclic RW1 for seasonal components (s_t)
seasonalComp2<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}

#Intrinsic GMRF density for spatial components (u_i)
logIGMRF1<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}


#' Build transition probability matrix
#'
#' @param G12 Probability of jumping to state 2 if currently at state 1.
#' @param G21 probability of jumping to state 1 if currently in state 2.
#'
#' @return A transition probability matrix
#' @export
#'
#' @examples G(0.1, 0.2)
#'
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

#Crude estimates
crudeEst<- function(y, e_it){
  ydot<- colSums(y, na.rm = T)
  edot<- colSums(e_it, na.rm = T)
  logydot<- log(ydot)
  logedot<- log(edot)
  lambdadot<- logydot-logedot
  nlambdadot<- lambdadot[lambdadot != -Inf]
  x<- 1:ncol(y)
  lambdadot<- ifelse(lambdadot== -Inf, mean(nlambdadot), lambdadot)
  success <- tryCatch({
    loess_fit <- loess(lambdadot ~ x, span = 0.3)
    TRUE
  }, error = function(e) {
    FALSE
  })

  if(success){
    loess_fit <- loess(lambdadot ~ x, span = 0.3)
    smoothed <- predict(loess_fit)
    crudeS<- lambdadot - smoothed
    crudeR<- smoothed
    #crudeU<- log(rowSums(y/e_it)/sum(exp(crudeR+crudeS)))-mean(log(rowSums(y/e_it)/sum(exp(crudeR+crudeS))))
    crudeU<- rep(0, nrow(y))
  }else{
    crudeR<- rep(mean(lambdadot), ncol(y))
    crudeS<- lambdadot-mean(lambdadot)
    crudeU<- rep(0, nrow(y))
    }
  return(list(crudeR, crudeS, crudeU))
}

#Check B for figures
checkB<- function(df){
  success <- tryCatch({
    Bs <- df$draws(variables = "B")
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(success)
}

#Extract posterior credible interval
posterior_interval_custom <- function(posterior_samples, prob = 0.95) {
  # Check if the input is a matrix
  if (!is.matrix(posterior_samples)) {
    stop("Input must be a matrix where rows are samples and columns are variables.")
  }

  lower <- (1 - prob) / 2
  upper <- 1 - lower

  # Apply quantile function across columns
  credible_intervals <- apply(posterior_samples, 2, quantile, probs = c(lower, upper))

  # Convert to dataframe and label rows
  credible_intervals_df <- as.data.frame(t(credible_intervals))
  colnames(credible_intervals_df) <- c("2.5%", "97.5%")

  return(credible_intervals_df)
}


#custom legend
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#Infer missing data for model evidence computation
InferredData<- function(list.ye, inf.object, adjmat, Model, burn.in = 2000){
  y<- list.ye[[1]]
  e_it<- list.ye[[2]]
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
  pred_Y<- matrix(NA, nrow = ndept, ncol = time)

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

  Alldata<- matrix(0, nrow = ndept, ncol = time)

  if(Model == 0){
    for(index in 1:length(thinning)){
      for(i in 1:ndept){
        for(t in 1:time){
          m<- (t - 1) %% 12 + 1
          Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
          pred_Y[i, t]<- rpois(1, Exlambda_it)
        }
      }
      Alldata<- Alldata + pred_Y
    }
  }else if(Model %in% c(1,2,4,5,7)){
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
      Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B.draws[index], Model = Model, adjmat = adjmat, Cyclic = T)
      for(i in 1:ndept){
        for(t in 1:time){
          m<- (t - 1) %% 12 + 1
          P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
          Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
          pred_Y[i, t]<- rpois(1, Exlambda_it)
        }
      }
      Alldata<- Alldata + pred_Y
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
      Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B, Model = Model, adjmat = adjmat, Cyclic = T)
      sum_Xit<- sum_Xit + Ex_Xit
      for(i in 1:ndept){
        for(t in 1:time){
          m<- (t - 1) %% 12 + 1
          P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
          Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index, 1] + P_Xit * z_it2[i, t] * B.draws[index, 2])
          pred_Y[i, t]<- rpois(1, Exlambda_it)
        }
      }
      Alldata<- Alldata + pred_Y
    }
  }
  Alldata<- round(Alldata/length(thinning))
  return(Alldata)
}

#Input missing data for model evidence computation
Inputdata<- function(y, e_it, inf.object, adjmat, Model){
  if(any(is.na(y))){
    ye<- list(y, e_it)
    inf.y<- InferredData(ye,inf.object,adjmat,Model)
    y<- ifelse(is.na(y), inf.y, y)
  }else{
    y<- y
  }
  return(y)
}

#Linear interpolation for application
LinearInterp<- function(popdata){
  ndept<- nrow(popdata)
  time<- ncol(popdata)*12
  e_it<- matrix(NA, nrow = ndept, ncol = time)
  availabledat<- seq(from = 1, to = time, by = 12)
  e_it[ , availabledat]<- popdata[ , 1:ncol(popdata)]

  for(i in 1:ndept){
    for(t in 1:(ncol(popdata)-1)){
      ind<- availabledat[t]
      ind2<- availabledat[t+1]
      df<- data.frame(x = c(ind,ind2), y = e_it[i, c(ind, ind2)])
      interv<- seq(from = ind+1, to = ind2-1, by = 1)
      e_it[i, (ind+1):(ind2-1)]<- approx(df$x, df$y, xout = interv)$y
    }
  }
  e_it<- e_it[ , -(max(availabledat):time)]
  return(round(e_it))
}


#Datasets for application (Meningococcal)
#Susceptibles<- read.csv("popn.csv")
#countries<- Susceptibles[ , 1]
#eit<- Susceptibles
#names(eit)<- NULL
#eit<- eit[ , -1]
#eit<- as.matrix(eit)
#eit<- LinearInterp(eit)
#dim(eit)
#Infected<- read.csv("Cases.csv")
#y<- matrix(Infected[ ,3], nrow = 28, ncol = 252, byrow = TRUE)
#alldata<- list(y, eit, countries)
#sim.plot(alldata)

#AdjacencyMatrix for application
#poly <- cshapes::cshp(date=as.Date("2019-12-31"), useGW=TRUE)
#Allcountriesnamescodes<- data.frame(countryname=poly$country_name, countrycode=poly$gwcode)
#requiredcountriesnames<-c("Austria","Belgium","Cyprus","Czech Republic",
#                         "Denmark","Estonia","Finland","France","German Federal Republic",
#                         "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
#                         "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
#                         "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom")

#1000km for 30%, 820km for 20%
#dmat <- cshapes::distmatrix(as.Date("2019-12-31"), type="capdist")
#colnames(dmat)<- Allcountriesnamescodes$countryname
#rownames(dmat)<- Allcountriesnamescodes$countryname
#dmat<- dmat[requiredcountriesnames, requiredcountriesnames]
#AdjacencyMatrix<- ifelse(dmat>820,0,1)
#diag(AdjacencyMatrix)<- 0
