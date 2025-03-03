
#' A function to compute the marginal log-likelihood for model comparison purposes in the spatio-temporal models presented in (Adeoye, et.al., 2025).
#'
#' @param y A space-time data matrix (locations on the row, and time on the columns).
#' @param e_it A space-time data matrix showing the number of susceptible individuals.
#' @param adjmat The adjacency matrix describing the connectivity of spatial locations.
#' @param Model The model specification (ranges from 0 to 6).
#' @param inf.object Resulting fit object from "infer2".
#' @param num_samples Number of samples to be used.
#'
#' @return The marginal log-likelihood
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' fit<-infer(TRUTHS[[1]],TRUTHS[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=FALSE,OutbreakProb=FALSE)
#' ModelEvidence(y=TRUTHS[[1]],e_it=TRUTHS[[2]],Model=0,adjmat=sim_adjmat,inf.object=fit)
#'
ModelEvidence<- function(y, e_it, adjmat, Model, inf.object, num_samples = 50000){

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  y<- ifelse(is.na(y), -1, y)

  time<- ncol(y)
  ndept<- nrow(y)

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  #R MCMC samples
  if(is.data.frame(inf.object)){
    if(Model == 0){
      inf.object<- inf.object[,-(c(1,2,(ncol(inf.object)-2):ncol(inf.object)))]
    }else if(Model %in% c(1, 2, 4, 5, 7)){
      inf.object<- inf.object[,-((ncol(inf.object)-1):ncol(inf.object))]
    }else if(Model %in% c(3, 6)){
      inf.object<- inf.object[, -ncol(inf.object)]
    }
  }else{
    #Stan HMC samples
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    if(Model == 0){
      inf.object<- inf.object[,-(1:2)]
    }
  }

  mu<- colMeans(inf.object)
  varcov<- cov(inf.object)
  cond<- check_cholesky(varcov)
  a<- numeric(num_samples)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
  }
  validDensities <- a[a != -Inf]
  MarginalLikelihood<- log(1) - log(length(validDensities)) + logSumExp_cpp(validDensities)
  return(MarginalLikelihood)
}


#' A function to compute the marginal log-likelihood for model comparison purposes in the spatio-temporal models presented in (Knorr-Held, et.al., 2003).
#'
#' @param y A space-time data matrix (locations on the row, and time on the columns).
#' @param e_it A space-time data matrix showing the number of susceptible individuals.
#' @param adjmat The adjacency matrix describing the connectivity of spatial locations.
#' @param Model The model specification (ranges from 0 to 6).
#' @param inf.object Resulting fit object from "infer".
#' @param num_samples Number of samples to be used.
#'
#' @return The marginal log-likelihood
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS1<- simulate.old(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' fit1<-infer.old(TRUTHS1[[1]],TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=FALSE,OutbreakProb=FALSE)
#' ModelEvidence.old(y=TRUTHS1[[1]],e_it=TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,inf.object=fit1)
#'
ModelEvidence.old<- function(y, e_it, adjmat, Model, inf.object, num_samples = 10000){

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  time<- ncol(y)
  ndept<- nrow(y)

  #R MCMC samples
  if(is.data.frame(inf.object)){
    if(Model == 0){
      inf.object<- inf.object[,-(c(1,2,(ncol(inf.object)-2):ncol(inf.object)))]
    }else if(Model %in% c(1, 2, 4, 5, 7)){
      inf.object<- inf.object[,-((ncol(inf.object)-1):ncol(inf.object))]
    }else if(Model %in% c(3, 6)){
      inf.object<- inf.object[, -ncol(inf.object)]
    }
  }else{
    #Stan HMC samples
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    if(Model == 0){
      inf.object<- inf.object[,-(1:2)]
    }
  }

  mu<- colMeans(inf.object)
  varcov<- cov(inf.object)

  cond<- check_cholesky(varcov)
  a<- numeric(num_samples)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:time)]
          u<- theta[3+time+time+(1:ndept)]
          ARcoeff<- 0
          a[i]<- GeneralLoglikelihood_cpp(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp(s, kappaS) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:time)]
          u<- theta[5+time+time+(1:ndept)]
          if(Model %in% c(3,6)){
            ARcoeff<- theta[5+time+time+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+time+ndept+1]
          }
          a[i]<- GeneralLoglikelihood_cpp(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp(s, kappaS) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:time)]
          u<- theta[3+time+time+(1:ndept)]
          ARcoeff<- 0
          a[i]<- GeneralLoglikelihood_cpp(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp(s, kappaS) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:time)]
          u<- theta[5+time+time+(1:ndept)]
          if(Model %in% c(3,6)){
            ARcoeff<- theta[5+time+time+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+time+ndept+1]
          }
          a[i]<- GeneralLoglikelihood_cpp(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp(s, kappaS) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
  }
  validDensities <- a[a != -Inf]
  MarginalLikelihood<- log(1) - log(length(validDensities)) + logSumExp_cpp(validDensities)
  return(MarginalLikelihood)
}
