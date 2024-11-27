#' A function to carryout Bayesian inference for the spatio-temporal models presented in (Adeoye, et.al., 2024).
#'
#' @param y A space-time data matrix (locations on the row, and time on the columns).
#' @param e_it A space-time data matrix showing the number of susceptible individuals (locations on the row, and time on the columns).
#' @param Model The model specification (ranges from 0 to 7).
#' @param adjmat The adjacency matrix describing the connectivity/neighborhood structure of the spatial locations.
#' @param num_iteration The required number of iterations for inference via MCMC. Default set to 30000.
#' @param Stan A logical argument asking whether to fit the model using Stan's HMC via cmdstanr.
#' @param nchains The number of parallel MCMC chains to run using Stan's HMC.
#' @param iter The required number of iterations for inference via HMC. Default set to 4000.
#' @param seed A random number for reproducibility.
#' @param verbose A logical argument asking whether to display (unnecessary) outputs.
#' @param ModEvid A logical argument asking whether to compute the model evidence.
#' @param OutbreakProb A logical argument asking whether to compute and return the posterior marginal probabilities of outbreak.
#' @param Burn.in MCMC burn-in period. Default set to 3000.
#' @param GPU A logical argument for accelerating Stan's computations using OpenCL (works with GPU only).
#' @param adaptdelta Target average proposal acceptance probability during Stan's adaptation period. Default set to 0.90.
#'
#' @return Samples from the target posterior distribution (and outbreak probabilities, model evidence).
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' HMCfit<- infer(TRUTHS[[1]],TRUTHS[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=TRUE)
#' HMCfit<- HMCfit[[1]]
#' mcmc.plot(HMCfit)
#' inf.plot(HMCfit)
#'
#' MCMCfit<- infer(TRUTHS[[1]],TRUTHS[[2]],Model=0,adjmat=sim_adjmat,Stan=FALSE,ModEvid=TRUE,OutbreakProb=TRUE)
#' MCMCfit<- MCMCfit[[1]]
#' mcmc.plot(MCMCfit)
#' inf.plot(MCMCfit)
#'
infer<- function(y, e_it, Model, adjmat, num_iteration = 30000, Stan = TRUE, GPU = FALSE, nchains = 4,
                  iter = 4000, seed = NULL, verbose = F, ModEvid = F, OutbreakProb = F, adaptdelta = 0.90, Burn.in = 1000){
  start_time <- Sys.time()

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  nstate<- 2
  ndept<- nrow(y)
  time<- ncol(y)

  original.y<- y

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }
  strs<- RW1PrecMat

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(original.y, adjmat)[[1]]
  z_it2 <- design_matrix_func(original.y, adjmat)[[2]]

  if(Model == 0){
    npar <- 0
  }
  else if(Model %in% c(1,2,4,5,7)){
    npar <- 1
  }
  else if(Model %in% c(3,6)){
    npar <- 2
  }

  if(Stan){
    initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), r = rep(0, time), sraw = rep(0, 11), kappa_u=20, kappa_r=20, kappa_s=20)
    initials_list <- lapply(1:nchains, function(x) initials)
    if(GPU){
      mod <- cmdstanr::cmdstan_model(system.file("stan", "newModel.stan", package = "DetectOutbreaks", mustWork = T), compile = F, cpp_options = list(stan_opencl = TRUE))
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                     SMat = strs, Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta, opencl_ids = c(0, 0))
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                       SMat = strs, Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta, opencl_ids = c(0, 0))
        })))
      }
    }else{
      mod <- cmdstanr::cmdstan_model(system.file("stan", "newModel.stan", package = "DetectOutbreaks", mustWork = T), compile = F)
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                     SMat = strs, Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta)
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                       SMat = strs, Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta)
        })))
      }
    }
    if(ModEvid){
      ME<- ModelEvidence(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = fit)
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = fit, adjmat = adjmat, Model = Model, Cyclic = T)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(fit, OutP))
    }else{
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(fit)
    }
  }else{

    crudeResults<- crudeEst(y, e_it)
    crudeR<- crudeResults[[1]]
    crudeS<- crudeResults[[2]]
    crudeU<- crudeResults[[3]]

    MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+2+1)
    initG12<- runif(1)
    initG21<- runif(1)
    initstateD<- state_dist_cpp(initG12, initG21)[2]
    MC_chain[1,]<- c(initG12, initG21, runif(1, 0, 1000), runif(1, 0, 1000), runif(1, 0, 30), crudeR, crudeS[1:12], crudeU, rep(0, 2), initstateD)

    zigmaR<- diag(rep(0.1, time), nrow = time, ncol = time)
    zigmaS<- diag(rep(0.1, 11), nrow = 11, ncol = 11)
    zigmaU<- diag(rep(0.08, ndept-1), nrow=ndept-1, ncol=ndept-1)
    optconstantR<- 2.38^2/(time-2)
    optconstantS<- 2.38^2/11
    optconstantU<- 2.38^2/(ndept-1)
    lambdaR<- 1
    lambdaS<- 1
    lambdaU<- 1

    RW2PrecMat<- matrix(0, nrow=time, ncol=time)
    RW2PrecMat[1,(1:3)]<- c(1,-2,1)
    RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
    RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
    RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
    RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
    for(i in 3:(time-3)){
      RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
    }
    strr<- RW2PrecMat

    Rsize<- floor(time/8)

    A<- 5+(1:Rsize)
    B<- 5+((Rsize+1):(2*Rsize))
    C<- 5+((2*Rsize+1):(3*Rsize))
    D<- 5+((3*Rsize+1):(4*Rsize))
    E<- 5+((4*Rsize+1):(5*Rsize))
    f<- 5+((5*Rsize+1):(6*Rsize))
    g<- 5+((6*Rsize+1):(7*Rsize))
    H<- 5+((7*Rsize+1):time)

    Blocks<- list(A,B,C,D,E,f,g,H)

    for(i in 2:num_iteration){
      #print(i)

      proposedkappaR<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% strr %*% MC_chain[i-1, 5+(1:time)])/2)
      MC_chain[i,3]<- proposedkappaR
      #print(paste("GibbskappaR = ", proposedkappaR))

      proposedkappaS<- rgamma(1, shape = 1 + 11/2, rate = (time / 12.0) * 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% strs %*% MC_chain[i-1, 5+time+(1:12)])/2)
      MC_chain[i,4]<- proposedkappaS
      #print(paste("GibbskappaS = ", proposedkappaS))

      proposedkappaU<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)
      MC_chain[i,5]<- proposedkappaU
      #print(paste("GibbskappaU = ", proposedkappaU))

      proposedspatcomps<- mvnfast::rmvn(1, mu=MC_chain[i-1, 5+time+12+(1:(ndept-1))], sigma = zigmaU)
      proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))

      priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R)
      priorproposedUcomps<- logIGMRF1(proposedspatcomps, MC_chain[i, 5], R)

      proposalproposedcompsU<- mvnfast::dmvn(proposedspatcomps[-ndept], mu = MC_chain[i-1, 5+time+12+(1:(ndept-1))], sigma = zigmaU, log = TRUE)
      proposalcurrentcompsU<- mvnfast::dmvn(MC_chain[i-1, 5+time+12+(1:(ndept-1))], mu = proposedspatcomps[-ndept], sigma = zigmaU, log = TRUE)

      likelihoodcurrent<- GeneralLoglikelihood_cpp2(y,MC_chain[i-1,5+(1:time)],MC_chain[i-1,5+time+(1:12)],MC_chain[i-1,5+time+12+(1:ndept)],G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp2(y,MC_chain[i-1,5+(1:time)],MC_chain[i-1,5+time+(1:12)],proposedspatcomps,G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

      mh.ratioU<- exp(likelihoodproposed + priorproposedUcomps + proposalcurrentcompsU
                      - likelihoodcurrent - priorcurrentUcomps - proposalproposedcompsU)

      #print(mh.ratioU)

      if(!is.na(mh.ratioU) && runif(1) < mh.ratioU){
        MC_chain[i,5+time+12+(1:ndept)]<- proposedspatcomps
      }
      else{
        MC_chain[i,5+time+12+(1:ndept)]<- MC_chain[i-1,5+time+12+(1:ndept)]
      }

      RW2PrecMat<- MC_chain[i, 3] * strr
      RconditionalcovA<- solve(RW2PrecMat[A-5, A-5])
      RconditionalcovB<- solve(RW2PrecMat[B-5, B-5])
      RconditionalcovC<- solve(RW2PrecMat[C-5, C-5])
      RconditionalcovD<- solve(RW2PrecMat[D-5, D-5])
      RconditionalcovE<- solve(RW2PrecMat[E-5, E-5])
      Rconditionalcovf<- solve(RW2PrecMat[f-5, f-5])
      Rconditionalcovg<- solve(RW2PrecMat[g-5, g-5])
      RconditionalcovH<- solve(RW2PrecMat[H-5, H-5])

      covBlocks<- list(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD,
                       RconditionalcovE, Rconditionalcovf, Rconditionalcovg, RconditionalcovH)

      for(j in 1:length(Blocks)){
        if(j==1){
          Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i-1, unlist(Blocks[-j])]
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(proposedRcomps, MC_chain[i-1, unlist(Blocks[-j])])
        }
        else if(j!=1 && j!=length(Blocks)){
          Rconditionalmean<- -covBlocks[[j]] %*% (RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[1:(j-1)]))-5] %*% MC_chain[i, unlist(Blocks[1:(j-1)])] + RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[(j+1):(length(Blocks))]))-5] %*% MC_chain[i, unlist(Blocks[(j+1):length(Blocks)])])
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(MC_chain[i, unlist(Blocks[1:(j-1)])], proposedRcomps, MC_chain[i, unlist(Blocks[(j+1):length(Blocks)])])
        }
        else if(j==length(Blocks)){
          Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i, unlist(Blocks[-j])]
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(MC_chain[i, unlist(Blocks[-j])], proposedRcomps)
        }
        likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i-1, 5+(1:time)], MC_chain[i-1, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
        likelihoodproposed<- GeneralLoglikelihood_cpp2(y, proposedRcomps, MC_chain[i-1, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)

        mh.ratioR<- exp(likelihoodproposed - likelihoodcurrent)

        #print(paste("mh.ratioR condpriorprop =", mh.ratioR))
        if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
          MC_chain[i, 5+(1:time)]<- proposedRcomps
        }
        else{
          if(j==1){
            MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
          }
          else if(j!=1){
            MC_chain[i, 5+(1:time)]<- MC_chain[i, 5+(1:time)]
          }
        }
      }

      proposedScomps<- mvnfast::rmvn(1, mu = MC_chain[i-1, 5+time+(1:11)], sigma = zigmaS)
      proposedScomps<- c(proposedScomps, -sum(proposedScomps))

      priorcurrentScomps<- seasonalComp2(MC_chain[i-1, 5+time+(1:12)], MC_chain[i, 4], strs)
      priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], strs)

      proposalproposedScomps<- mvnfast::dmvn(proposedScomps[-12], mu = MC_chain[i-1, 5+time+(1:11)], sigma = zigmaS, log = TRUE)
      proposalcurrentScomps<- mvnfast::dmvn(MC_chain[i-1, 5+time+(1:11)], mu = proposedScomps[-12], sigma = zigmaS, log = TRUE)

      likelihoodcurrent<- GeneralLoglikelihood_cpp2(y,MC_chain[i,5+(1:time)],MC_chain[i-1,5+time+(1:12)],MC_chain[i,5+time+12+(1:ndept)],G(MC_chain[i-1,1], MC_chain[i-1,2]),e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp2(y,MC_chain[i,5+(1:time)],proposedScomps,MC_chain[i,5+time+12+(1:ndept)],G(MC_chain[i-1,1], MC_chain[i-1,2]),e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

      mh.ratioS<- exp(likelihoodproposed + priorproposedScomps + proposalcurrentScomps
                      - likelihoodcurrent - priorcurrentScomps - proposalproposedScomps)

      #print(mh.ratioS)

      if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
        MC_chain[i,5+time+(1:12)]<- proposedScomps
      }
      else{
        MC_chain[i,5+time+(1:12)]<- MC_chain[i-1,5+time+(1:12)]
      }

      if(Model == 0){
        proposedB <- c(0, 0)
        priorcurrentB<- -99999
        priorproposedB<- -99999
      }else if(Model %in% c(1,2,4,5,7)) {
        proposedB <- abs(rnorm(1, mean = MC_chain[i-1, 5+time+12+ndept+1], sd = 0.1))
        proposedB <- c(proposedB, 0)
        priorcurrentB<- dgamma(MC_chain[i-1, 5+time+12+ndept+1], shape = 2, rate = 2, log=TRUE)
        priorproposedB<- dgamma(proposedB[1], shape = 2, rate = 2, log=TRUE)
      }else if(Model %in% c(3,6)){
        proposedB <- abs(rnorm(2, mean = MC_chain[i-1, 5+time+12+ndept+(1:2)], sd = c(0.1, 0.1)))
        priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:2)], shape = c(2, 2), rate = c(2,2), log=TRUE))
        priorproposedB<- sum(dgamma(proposedB, shape = c(2, 2), rate = c(2, 2), log=TRUE))
      }

      likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1],MC_chain[i-1,2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1,2]), e_it, proposedB, Model,z_it, z_it2)

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 5+time+12+ndept+(1:2)]<- proposedB
      }
      else{
        MC_chain[i, 5+time+12+ndept+(1:2)]<- MC_chain[i-1, 5+time+12+ndept+(1:2)]
      }

      proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(0.1, 0.1)))
      if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
      if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
      priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

      likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i,5+(1:time)], MC_chain[i,5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(proposedGs[1], proposedGs[2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

      mh.ratio<- exp(likelihoodproposed + priorproposedGs
                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 1:2]<- proposedGs
      }
      else{
        MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
      }
      MC_chain[i, 5+time+12+ndept+2+1]<- state_dist_cpp(MC_chain[i, 1], MC_chain[i, 2])[2]

      #Adapting zigmaR
      if(i==5){
        epsilonR<- 0.007
        XnR<- MC_chain[1:i, 5+(1:time)]
        XnbarR <- colMeans(XnR)
        zigmaR <- cov(XnR) + epsilonR * diag(rep(1, time))
        zigmaR<- optconstantR * zigmaR
      } else if (i > 5){

        ### Using random walk after 5 conditional prior proposals
        proposedRcomps<- mvnfast::rmvn(1, mu = MC_chain[i, 5+(1:time)], sigma = zigmaR)

        priorcurrentRcomps<- randomwalk2(MC_chain[i, 5+(1:time)], MC_chain[i, 3])
        priorproposedRcomps<- randomwalk2(proposedRcomps, MC_chain[i, 3])

        proposalproposedRcomps<- mvnfast::dmvn(proposedRcomps, mu = MC_chain[i, 5+(1:time)], sigma = zigmaR, log = TRUE)
        proposalcurrentRcomps<- mvnfast::dmvn(MC_chain[i, 5+(1:time)], mu = proposedRcomps, sigma = zigmaR, log = TRUE)

        likelihoodcurrent<- GeneralLoglikelihood_cpp2(y,MC_chain[i,5+(1:time)],MC_chain[i, 5+time+(1:12)],MC_chain[i,5+time+12+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)
        likelihoodproposed<- GeneralLoglikelihood_cpp2(y,proposedRcomps,MC_chain[i, 5+time+(1:12)],MC_chain[i,5+time+12+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

        mh.ratioR<- exp(likelihoodproposed + priorproposedRcomps + proposalcurrentRcomps
                        - likelihoodcurrent - priorcurrentRcomps - proposalproposedRcomps)

        #print(mh.ratioR)

        if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
          MC_chain[i,5+(1:time)]<- proposedRcomps
        }
        else{
          MC_chain[i,5+(1:time)]<- MC_chain[i,5+(1:time)]
        }

        XnbarPrevR <- XnbarR
        XnbarR <- (i*XnbarR + MC_chain[i, 5+(1:time)])/(i+1)
        zigmaR <- ((i-1)*zigmaR + tcrossprod(MC_chain[i, 5+(1:time)]) + i*tcrossprod(XnbarPrevR) - (i+1)*tcrossprod(XnbarR) + epsilonR*diag(rep(1,time)))/i
        #Robbins Munro tuning
        lambdaR<- lambdaR * exp((2/max(1, i-5)) * (min(mh.ratioR, 1) - 0.234))
        zigmaR<- lambdaR* optconstantR * zigmaR
        #print(zigmaR)
      }

      #Adapting zigmaS
      if(i==5){
        epsilonS<- 0.007
        XnS<- MC_chain[1:i, 5+time+(1:11)]
        XnbarS <- colMeans(XnS)
        zigmaS <- cov(XnS) + epsilonS*diag(rep(1, 11))
        zigmaS<- optconstantS * zigmaS
      } else if (i > 5){
        XnbarPrevS <- XnbarS
        XnbarS <- (i*XnbarS + MC_chain[i, 5+time+(1:11)])/(i+1)
        zigmaS <- ((i-1)*zigmaS + tcrossprod(MC_chain[i, 5+time+(1:11)]) + i*tcrossprod(XnbarPrevS) - (i+1)*tcrossprod(XnbarS) + epsilonS*diag(rep(1,11)))/i
        #Robbins Munro tuning
        lambdaS<- lambdaS * exp((2/max(1, i-5)) * (min(mh.ratioS, 1) - 0.234))
        zigmaS<- lambdaS* optconstantS * zigmaS
        #print(zigmaS)
      }

      #Adapting zigmaU
      if(i==5){
        epsilonU<- 0.007
        XnU<- MC_chain[1:i, 5+time+12+(1:(ndept-1))]
        XnbarU <- colMeans(XnU)
        zigmaU <- cov(XnU) + epsilonU*diag(rep(1, ndept-1))
        zigmaU<- optconstantU * zigmaU
      } else if (i > 5){
        XnbarPrevU <- XnbarU
        XnbarU <- (i*XnbarU + MC_chain[i, 5+time+12+(1:(ndept-1))])/(i+1)
        zigmaU <- ((i-1)*zigmaU + tcrossprod(MC_chain[i, 5+time+12+(1:(ndept-1))]) + i*tcrossprod(XnbarPrevU) - (i+1)*tcrossprod(XnbarU) + epsilonU*diag(rep(1,ndept-1)))/i
        #Robbins Munro tuning
        lambdaU<- lambdaU * exp((2/max(1, i-5)) * (min(mh.ratioU, 1) - 0.234))
        zigmaU<- lambdaU* optconstantU * zigmaU
        #print(zigmaU)
      }
    }

    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:2, sep=""), "StationaryDistribution"))
    MC_chain<- as.data.frame(MC_chain)
    if(ModEvid){
      ME<- ModelEvidence(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = MC_chain[-(1:2000), ])
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = MC_chain, adjmat = adjmat, Model = Model, Cyclic = T)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(MC_chain, OutP))
    }else{
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(MC_chain)
    }
  }
}



#' A function to carryout Bayesian inference for the spatio-temporal models presented in (Knorr-Held, et.al., 2003).
#'
#' @param y A space-time data matrix (locations on the rows, and time on the columns).
#' @param e_it A space-time data matrix showing the number of susceptible individuals (locations on the rows, and time on the columns).
#' @param Model The model specification (ranges from 0 to 7).
#' @param adjmat The adjacency matrix describing the connectivity/neighborhood structure of the spatial locations.
#' @param num_iteration The required number of iterations for inference via MCMC. Default set to 30000.
#' @param Stan A logical argument asking whether to fit the model using Stan's HMC via cmdstanr.
#' @param nchains The number of parallel MCMC chains to run using Stan's HMC.
#' @param iter The required number of iterations for inference via HMC.
#' @param seed A random number for reproducibility.
#' @param verbose A logical argument asking whether to display (unnecessary) outputs.
#' @param ModEvid A logical argument asking whether to compute and print the model evidence.
#' @param OutbreakProb A logical argument asking whether to compute and return the posterior marginal probabilities of outbreak.
#' @param Burn.in MCMC burn-in period. Default set to 3000.
#' @param GPU A logical argument for accelerating Stan's computations using OpenCL (works with GPU only).
#' @param adaptdelta Target average proposal acceptance probability during Stan's adaptation period. Default set to 0.90.
#'
#' @return Samples from the target posterior distribution (and outbreak probabilities, model evidence).
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS1<- simulate.old(Model = 0, time = 48, adj.matrix = sim_adjmat);
#' HMCfit1<-infer.old(TRUTHS1[[1]],TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=TRUE);
#' HMCfit1<- HMCfit1[[1]]
#' mcmc.plot(HMCfit1)
#' inf.plot(HMCfit1)
#'
#'MCMCfit1<-infer.old(TRUTHS1[[1]],TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,Stan=FALSE,ModEvid=TRUE,OutbreakProb=TRUE);
#'
infer.old<- function(y, e_it, Model, adjmat, num_iteration = 30000, Stan = TRUE, GPU = FALSE, nchains = 4,
                  iter = 4000, seed = NULL, verbose = F, ModEvid = F, OutbreakProb = F, adaptdelta = 0.90, Burn.in = 1000){
  start_time <- Sys.time()

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  nstate<- 2
  ndept<- nrow(y)
  time<- ncol(y)

  original.y<- y

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(original.y, adjmat)[[1]]
  z_it2 <- design_matrix_func(original.y, adjmat)[[2]]

  if(Model == 0){
    npar <- 0
  }
  else if(Model %in% c(1,2,4,5,7)){
    npar <- 1
  }
  else if(Model %in% c(3,6)){
    npar <- 2
  }

  if(Stan){
    initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), r = rep(0, time), s = rep(0, time), kappa_u=20, kappa_r=20, kappa_s=20)
    initials_list <- lapply(1:nchains, function(x) initials)
    if(GPU){
      mod <- cmdstanr::cmdstan_model(system.file("stan", "Model.stan", package = "DetectOutbreaks", mustWork = T), compile = F, cpp_options = list(stan_opencl = TRUE))
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                     Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta, opencl_ids = c(0, 0))
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                       Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta, opencl_ids = c(0, 0))
        })))
      }
    }else{
      mod <- cmdstanr::cmdstan_model(system.file("stan", "Model.stan", package = "DetectOutbreaks", mustWork = T), compile = F)
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                     Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta)
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, y=y, e_it=e_it, R=R,
                                       Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta)
        })))
      }
    }
    if(ModEvid){
      ME<- ModelEvidence.old(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = fit)
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = fit, adjmat = adjmat, Model = Model, Cyclic = F)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(fit, OutP))
    }else{
    end_time <- Sys.time()
    time_taken<- end_time - start_time
    print(time_taken)
    return(fit)
    }
  }else{

    crudeResults<- crudeEst(y, e_it)
    crudeR<- crudeResults[[1]]
    crudeS<- crudeResults[[2]]
    crudeU<- crudeResults[[3]]

    MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+time+ndept+2+1)
    initG12<- runif(1)
    initG21<- runif(1)
    initstateD<- state_dist_cpp(initG12, initG21)[2]
    MC_chain[1,]<- c(initG12, initG21, runif(1, 0, 1000), runif(1, 0, 1000), runif(1, 0, 30), crudeR, crudeS, crudeU, rep(0, 2), initstateD)

    zigmaR<- diag(rep(0.1, time), nrow = time, ncol = time)
    zigmaS<- diag(rep(0.1, time), nrow = time, ncol = time)
    zigmaU<- diag(rep(0.08, ndept-1), nrow=ndept-1, ncol=ndept-1)
    optconstantR<- 2.38^2/(time-2)
    optconstantS<- 2.38^2/(time-11)
    optconstantU<- 2.38^2/(ndept-1)
    lambdaR<- 1
    lambdaS<- 1
    lambdaU<- 1

    RW2PrecMat<- matrix(0, nrow=time, ncol=time)
    RW2PrecMat[1,(1:3)]<- c(1,-2,1)
    RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
    RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
    RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
    RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
    for(i in 3:(time-3)){
      RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
    }
    strr<- RW2PrecMat

    Rsize<- floor(time/8)

    A<- 5+(1:Rsize)
    B<- 5+((Rsize+1):(2*Rsize))
    C<- 5+((2*Rsize+1):(3*Rsize))
    D<- 5+((3*Rsize+1):(4*Rsize))
    E<- 5+((4*Rsize+1):(5*Rsize))
    f<- 5+((5*Rsize+1):(6*Rsize))
    g<- 5+((6*Rsize+1):(7*Rsize))
    H<- 5+((7*Rsize+1):time)

    Blocks<- list(A,B,C,D,E,f,g,H)

    SRWPrecMat<- matrix(0, nrow=time, ncol=time)
    SRWPrecMat[1,(1:12)]<- rep(1, 12)
    for(i in 2:12){
      SRWPrecMat[i, 1:(i+11)]<- c(1:(i-1), rep(i, (12+1-i)), (i-1):1)
    }
    for(i in 14:(time-11)){
      SRWPrecMat[i-1, ((i-12):(i+10))]<- c(1:12,11:1)
    }
    for(i in (time-1):(time-11)){
      SRWPrecMat[i, (i-11):time]<- c(1:(time-i), rep((time-i)+1, i-(time-12)), (time-i):1)
    }
    SRWPrecMat[time, ((time-11):time)]<- rep(1, 12)
    strs<- SRWPrecMat

    size<- floor(time/12)

    SA<- 5+time+(1:size)
    SB<- 5+time+((size+1):(2*size))
    SC<- 5+time+((2*size+1):(3*size))
    SD<- 5+time+((3*size+1):(4*size))
    SE<- 5+time+((4*size+1):(5*size))
    Sf<- 5+time+((5*size+1):(6*size))
    Sg<- 5+time+((6*size+1):(7*size))
    SH<- 5+time+((7*size+1):(8*size))
    SI<- 5+time+((8*size+1):(9*size))
    SJ<- 5+time+((9*size+1):(10*size))
    SK<- 5+time+((10*size+1):(11*size))
    SL<- 5+time+((11*size+1):time)

    SBlocks<- list(SA,SB,SC,SD,SE,Sf,Sg,SH,SI,SJ,SK,SL)

    for(i in 2:num_iteration){
      #print(i)

      proposedkappaR<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% strr %*% MC_chain[i-1, 5+(1:time)])/2)
      MC_chain[i,3]<- proposedkappaR
      #print(paste("GibbskappaR = ", proposedkappaR))

      proposedkappaS<- rgamma(1, shape = 1 + (time-11)/2, rate = 0.001 + (t(MC_chain[i-1, 5+time+(1:time)]) %*% strs %*% MC_chain[i-1, 5+time+(1:time)])/2)
      MC_chain[i,4]<- proposedkappaS
      #print(paste("GibbskappaS = ", proposedkappaS))

      proposedkappaU<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+time+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+time+(1:ndept)])/2)
      MC_chain[i,5]<- proposedkappaU
      #print(paste("GibbskappaU = ", proposedkappaU))

      proposedspatcomps<- mvnfast::rmvn(1, mu=MC_chain[i-1, 5+time+time+(1:(ndept-1))], sigma = zigmaU)
      proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))

      priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+time+(1:ndept)], MC_chain[i, 5], R)
      priorproposedUcomps<- logIGMRF1(proposedspatcomps, MC_chain[i, 5], R)

      proposalproposedcompsU<- mvnfast::dmvn(proposedspatcomps[-ndept], mu = MC_chain[i-1, 5+time+time+(1:(ndept-1))], sigma = zigmaU, log = TRUE)
      proposalcurrentcompsU<- mvnfast::dmvn(MC_chain[i-1, 5+time+time+(1:(ndept-1))], mu = proposedspatcomps[-ndept], sigma = zigmaU, log = TRUE)

      likelihoodcurrent<- GeneralLoglikelihood_cpp(y,MC_chain[i-1,5+(1:time)],MC_chain[i-1,5+time+(1:time)],MC_chain[i-1,5+time+time+(1:ndept)],G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp(y,MC_chain[i-1,5+(1:time)],MC_chain[i-1,5+time+(1:time)],proposedspatcomps,G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)

      mh.ratioU<- exp(likelihoodproposed + priorproposedUcomps + proposalcurrentcompsU
                      - likelihoodcurrent - priorcurrentUcomps - proposalproposedcompsU)

      #print(mh.ratioU)

      if(!is.na(mh.ratioU) && runif(1) < mh.ratioU){
        MC_chain[i,5+time+time+(1:ndept)]<- proposedspatcomps
      }
      else{
        MC_chain[i,5+time+time+(1:ndept)]<- MC_chain[i-1,5+time+time+(1:ndept)]
      }

      RW2PrecMat<- MC_chain[i, 3] * strr
      RconditionalcovA<- solve(RW2PrecMat[A-5, A-5])
      RconditionalcovB<- solve(RW2PrecMat[B-5, B-5])
      RconditionalcovC<- solve(RW2PrecMat[C-5, C-5])
      RconditionalcovD<- solve(RW2PrecMat[D-5, D-5])
      RconditionalcovE<- solve(RW2PrecMat[E-5, E-5])
      Rconditionalcovf<- solve(RW2PrecMat[f-5, f-5])
      Rconditionalcovg<- solve(RW2PrecMat[g-5, g-5])
      RconditionalcovH<- solve(RW2PrecMat[H-5, H-5])

      covBlocks<- list(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD,
                       RconditionalcovE, Rconditionalcovf, Rconditionalcovg, RconditionalcovH)

      SRWPrecMat<- MC_chain[i, 4] * strs
      SconditionalcovA<- solve(SRWPrecMat[SA-5-time, SA-5-time])
      SconditionalcovB<- solve(SRWPrecMat[SB-5-time, SB-5-time])
      SconditionalcovC<- solve(SRWPrecMat[SC-5-time, SC-5-time])
      SconditionalcovD<- solve(SRWPrecMat[SD-5-time, SD-5-time])
      SconditionalcovE<- solve(SRWPrecMat[SE-5-time, SE-5-time])
      Sconditionalcovf<- solve(SRWPrecMat[Sf-5-time, Sf-5-time])
      Sconditionalcovg<- solve(SRWPrecMat[Sg-5-time, Sg-5-time])
      SconditionalcovH<- solve(SRWPrecMat[SH-5-time, SH-5-time])
      SconditionalcovI<- solve(SRWPrecMat[SI-5-time, SI-5-time])
      SconditionalcovJ<- solve(SRWPrecMat[SJ-5-time, SJ-5-time])
      SconditionalcovK<- solve(SRWPrecMat[SK-5-time, SK-5-time])
      SconditionalcovL<- solve(SRWPrecMat[SL-5-time, SL-5-time])

      ScovBlocks<- list(SconditionalcovA, SconditionalcovB, SconditionalcovC, SconditionalcovD, SconditionalcovE, Sconditionalcovf, Sconditionalcovg,
                        SconditionalcovH, SconditionalcovI, SconditionalcovJ, SconditionalcovK, SconditionalcovL)

      for(j in 1:length(Blocks)){
        if(j==1){
          Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i-1, unlist(Blocks[-j])]
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(proposedRcomps, MC_chain[i-1, unlist(Blocks[-j])])
        }
        else if(j!=1 && j!=length(Blocks)){
          Rconditionalmean<- -covBlocks[[j]] %*% (RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[1:(j-1)]))-5] %*% MC_chain[i, unlist(Blocks[1:(j-1)])] + RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[(j+1):(length(Blocks))]))-5] %*% MC_chain[i, unlist(Blocks[(j+1):length(Blocks)])])
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(MC_chain[i, unlist(Blocks[1:(j-1)])], proposedRcomps, MC_chain[i, unlist(Blocks[(j+1):length(Blocks)])])
        }
        else if(j==length(Blocks)){
          Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i, unlist(Blocks[-j])]
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(MC_chain[i, unlist(Blocks[-j])], proposedRcomps)
        }

        likelihoodcurrent<- GeneralLoglikelihood_cpp(y, MC_chain[i-1, 5+(1:time)], MC_chain[i-1, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)
        likelihoodproposed<- GeneralLoglikelihood_cpp(y, proposedRcomps, MC_chain[i-1, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)

        mh.ratioR<- exp(likelihoodproposed - likelihoodcurrent)

        #print(paste("mh.ratioR condpriorprop =", mh.ratioR))
        if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
          MC_chain[i, 5+(1:time)]<- proposedRcomps
        }
        else{
          if(j==1){
            MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
          }
          else if(j!=1){
            MC_chain[i, 5+(1:time)]<- MC_chain[i, 5+(1:time)]
          }
        }
      }

      for(j in 1:length(SBlocks)){
        if(j==1){
          Sconditionalmean<- -ScovBlocks[[j]] %*% SRWPrecMat[(SBlocks[[j]])-(5+time), (unlist(SBlocks[-j]))-(5+time)] %*% MC_chain[i-1, unlist(SBlocks[-j])]
          proposedScomps<- mvnfast::rmvn(1, mu = Sconditionalmean, sigma = ScovBlocks[[j]])
          proposedScomps<- c(proposedScomps, MC_chain[i-1, unlist(SBlocks[-j])])
        }
        else if(j!=1 && j!=length(SBlocks)){
          Sconditionalmean<- -ScovBlocks[[j]] %*% (SRWPrecMat[(SBlocks[[j]])-(5+time), (unlist(SBlocks[1:(j-1)]))-(5+time)] %*% MC_chain[i, unlist(SBlocks[1:(j-1)])] + SRWPrecMat[(SBlocks[[j]])-(5+time), (unlist(SBlocks[(j+1):length(SBlocks)]))-(5+time)] %*% MC_chain[i, unlist(SBlocks[(j+1):length(SBlocks)])])
          proposedScomps<- mvnfast::rmvn(1, mu = Sconditionalmean, sigma = ScovBlocks[[j]])
          proposedScomps<- c(MC_chain[i, unlist(SBlocks[1:(j-1)])], proposedScomps, MC_chain[i, unlist(SBlocks[(j+1):length(SBlocks)])])
        }
        else if(j==length(SBlocks)){
          Sconditionalmean<- -ScovBlocks[[j]] %*% SRWPrecMat[(SBlocks[[j]])-(5+time), (unlist(SBlocks[-j]))-(5+time)] %*% MC_chain[i, unlist(SBlocks[-j])]
          proposedScomps<- mvnfast::rmvn(1, mu = Sconditionalmean, sigma = ScovBlocks[[j]])
          proposedScomps<- c(MC_chain[i, unlist(SBlocks[-j])], proposedScomps)
        }

        likelihoodcurrent<- GeneralLoglikelihood_cpp(y,MC_chain[i, 5+(1:time)],MC_chain[i-1,5+time+(1:time)],MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
        likelihoodproposed<- GeneralLoglikelihood_cpp(y,MC_chain[i, 5+(1:time)],proposedScomps,MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)

        mh.ratioS<- exp(likelihoodproposed - likelihoodcurrent)

        #print(paste("mh.ratioS condpriorprop = ", mh.ratioS))
        if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
          MC_chain[i, 5+time+(1:time)]<- proposedScomps
        }
        else{
          if(j==1){
            MC_chain[i, 5+time+(1:time)]<- MC_chain[i-1, 5+time+(1:time)]
          }
          else if(j!=1){
            MC_chain[i, 5+time+(1:time)]<- MC_chain[i, 5+time+(1:time)]
          }
        }
      }

      if(Model == 0){
        proposedB <- c(0, 0)
        priorcurrentB<- -99999
        priorproposedB<- -99999
      }else if(Model %in% c(1,2,4,5,7)) {
        proposedB <- abs(rnorm(1, mean = MC_chain[i-1, 5+time+time+ndept+1], sd = 0.1))
        proposedB <- c(proposedB, 0)
        priorcurrentB<- dgamma(MC_chain[i-1, 5+time+time+ndept+1], shape = 2, rate = 2, log=TRUE)
        priorproposedB<- dgamma(proposedB[1], shape = 2, rate = 2, log=TRUE)
      }else if(Model %in% c(3,6)){
        proposedB <- abs(rnorm(2, mean = MC_chain[i-1, 5+time+time+ndept+(1:2)], sd = c(0.1, 0.1)))
        priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+time+ndept+(1:2)], shape = c(2,2), rate = c(2,2), log=TRUE))
        priorproposedB<- sum(dgamma(proposedB, shape = c(2,2), rate = c(2, 2), log=TRUE))
      }

      likelihoodcurrent<- GeneralLoglikelihood_cpp(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1],MC_chain[i-1,2]), e_it, MC_chain[i-1, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1,2]), e_it, proposedB, Model,z_it, z_it2)

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 5+time+time+ndept+(1:2)]<- proposedB
      }
      else{
        MC_chain[i, 5+time+time+ndept+(1:2)]<- MC_chain[i-1, 5+time+time+ndept+(1:2)]
      }

      proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(0.1, 0.1)))
      if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
      if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
      priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

      likelihoodcurrent<- GeneralLoglikelihood_cpp(y, MC_chain[i,5+(1:time)], MC_chain[i,5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model, z_it, z_it2)
      likelihoodproposed<- GeneralLoglikelihood_cpp(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:time)], MC_chain[i, 5+time+time+(1:ndept)], G(proposedGs[1], proposedGs[2]),e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)

      mh.ratio<- exp(likelihoodproposed + priorproposedGs
                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 1:2]<- proposedGs
      }
      else{
        MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
      }
      MC_chain[i, 5+time+time+ndept+2+1]<- state_dist_cpp(MC_chain[i, 1], MC_chain[i, 2])[2]

      #Adapting zigmaR
      if(i==5){
        epsilonR<- 0.07
        XnR<- MC_chain[1:i, 5+(1:time)]
        XnbarR <- colMeans(XnR)
        zigmaR <- cov(XnR) + epsilonR * diag(rep(1, time))
        zigmaR<- optconstantR * zigmaR
      } else if (i > 5){

        ### Using random walk after 5 conditional prior proposals
        proposedRcomps<- mvnfast::rmvn(1, mu = MC_chain[i, 5+(1:time)], sigma = zigmaR)

        priorcurrentRcomps<- randomwalk2(MC_chain[i, 5+(1:time)], MC_chain[i, 3])
        priorproposedRcomps<- randomwalk2(proposedRcomps, MC_chain[i, 3])

        proposalproposedRcomps<- mvnfast::dmvn(proposedRcomps, mu = MC_chain[i, 5+(1:time)], sigma = zigmaR, log = TRUE)
        proposalcurrentRcomps<- mvnfast::dmvn(MC_chain[i, 5+(1:time)], mu = proposedRcomps, sigma = zigmaR, log = TRUE)

        likelihoodcurrent<- GeneralLoglikelihood_cpp(y,MC_chain[i,5+(1:time)],MC_chain[i, 5+time+(1:time)],MC_chain[i,5+time+time+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]),e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
        likelihoodproposed<- GeneralLoglikelihood_cpp(y,proposedRcomps,MC_chain[i, 5+time+(1:time)],MC_chain[i,5+time+time+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]),e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)

        mh.ratioR<- exp(likelihoodproposed + priorproposedRcomps + proposalcurrentRcomps
                        - likelihoodcurrent - priorcurrentRcomps - proposalproposedRcomps)

        #print(mh.ratioR)

        if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
          MC_chain[i,5+(1:time)]<- proposedRcomps
        }
        else{
          MC_chain[i,5+(1:time)]<- MC_chain[i,5+(1:time)]
        }

        proposedScomps<- mvnfast::rmvn(1, mu = MC_chain[i, 5+time+(1:time)], sigma = zigmaS)

        priorcurrentScomps<- seasonalComp(MC_chain[i, 5+time+(1:time)], MC_chain[i, 4])
        priorproposedScomps<- seasonalComp(proposedScomps, MC_chain[i, 4])

        proposalproposedScomps<- mvnfast::dmvn(proposedScomps, mu = MC_chain[i, 5+time+(1:time)], sigma = zigmaS, log = TRUE)
        proposalcurrentScomps<- mvnfast::dmvn(MC_chain[i, 5+time+(1:time)], mu = proposedScomps, sigma = zigmaS, log = TRUE)

        likelihoodcurrent<- GeneralLoglikelihood_cpp(y,MC_chain[i,5+(1:time)],MC_chain[i,5+time+(1:time)],MC_chain[i,5+time+time+(1:ndept)],G(MC_chain[i,1], MC_chain[i,2]),e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)
        likelihoodproposed<- GeneralLoglikelihood_cpp(y,MC_chain[i,5+(1:time)],proposedScomps,MC_chain[i,5+time+time+(1:ndept)],G(MC_chain[i,1], MC_chain[i,2]),e_it, MC_chain[i, 5+time+time+ndept+(1:2)], Model,z_it, z_it2)

        mh.ratioS<- exp(likelihoodproposed + priorproposedScomps + proposalcurrentScomps
                        - likelihoodcurrent - priorcurrentScomps - proposalproposedScomps)

        #print(mh.ratioS)

        if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
          MC_chain[i,5+time+(1:time)]<- proposedScomps
        }
        else{
          MC_chain[i,5+time+(1:time)]<- MC_chain[i,5+time+(1:time)]
        }

        XnbarPrevR <- XnbarR
        XnbarR <- (i*XnbarR + MC_chain[i, 5+(1:time)])/(i+1)
        zigmaR <- ((i-1)*zigmaR + tcrossprod(MC_chain[i, 5+(1:time)]) + i*tcrossprod(XnbarPrevR) - (i+1)*tcrossprod(XnbarR) + epsilonR*diag(rep(1,time)))/i
        #Robbins Munro tuning
        lambdaR<- lambdaR * exp((2/max(1, i-5)) * (min(mh.ratioR, 1) - 0.234))
        zigmaR<- lambdaR* optconstantR * zigmaR
        #print(zigmaR)
      }

      #Adapting zigmaS
      if(i==5){
        epsilonS<- 0.07
        XnS<- MC_chain[1:i, 5+time+(1:time)]
        XnbarS <- colMeans(XnS)
        zigmaS <- cov(XnS) + epsilonS*diag(rep(1, time))
        zigmaS<- optconstantS * zigmaS
      } else if (i > 5){
        XnbarPrevS <- XnbarS
        XnbarS <- (i*XnbarS + MC_chain[i, 5+time+(1:time)])/(i+1)
        zigmaS <- ((i-1)*zigmaS + tcrossprod(MC_chain[i, 5+time+(1:time)]) + i*tcrossprod(XnbarPrevS) - (i+1)*tcrossprod(XnbarS) + epsilonS*diag(rep(1,time)))/i
        #Robbins Munro tuning
        lambdaS<- lambdaS * exp((2/max(1, i-5)) * (min(mh.ratioS, 1) - 0.234))
        zigmaS<- lambdaS* optconstantS * zigmaS
        #print(zigmaS)
      }

      #Adapting zigmaU
      if(i==5){
        epsilonU<- 0.07
        XnU<- MC_chain[1:i, 5+time+time+(1:(ndept-1))]
        XnbarU <- colMeans(XnU)
        zigmaU <- cov(XnU) + epsilonU*diag(rep(1, ndept-1))
        zigmaU<- optconstantU * zigmaU
      } else if (i > 5){
        XnbarPrevU <- XnbarU
        XnbarU <- (i*XnbarU + MC_chain[i, 5+time+time+(1:(ndept-1))])/(i+1)
        zigmaU <- ((i-1)*zigmaU + tcrossprod(MC_chain[i, 5+time+time+(1:(ndept-1))]) + i*tcrossprod(XnbarPrevU) - (i+1)*tcrossprod(XnbarU) + epsilonU*diag(rep(1,ndept-1)))/i
        #Robbins Munro tuning
        lambdaU<- lambdaU * exp((2/max(1, i-5)) * (min(mh.ratioU, 1) - 0.234))
        zigmaU<- lambdaU* optconstantU * zigmaU
        #print(zigmaU)
      }
    }

    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:time, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:2, sep=""), "StationaryDistribution"))
    MC_chain<- as.data.frame(MC_chain)
    if(ModEvid){
      ME<- ModelEvidence(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = MC_chain[-(1:2000), ])
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = MC_chain, adjmat = adjmat, Model = Model, Cyclic = F)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(MC_chain, OutP))
    }else{
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(MC_chain)
    }
  }
}

