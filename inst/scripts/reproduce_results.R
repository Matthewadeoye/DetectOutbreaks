#Set up Stan
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
install_cmdstan()
# install our new package - DetectOutbreaks, and other packages for figures
devtools::install_github("Matthewadeoye/DetectOutbreaks")
packages<- c("gdata","pROC", "ggplot2", "sf", "dplyr", "RColorBrewer", "cshapes", "GGally")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
}))

#Set-up adjacency matrix for simulation-study (See Figure 1 connectivity in paper)
sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang

# Note: to see trace-plots for MCMC fits, simply do "DetectOutbreaks::mcmc.plot(fitobject)", where fitobject is the resulting HMC or MCMC fit.

#Section 3 - Application to simulated data
#Simulate from each model (See Table 1 and Table 2 in paper)
B<- list(c(0, 0), c(1.65, 0), c(1.65, 0), c(1.25, 0.75), c(0.55, 0), c(0.49, 0), c(0.35, 0.20), c(1.65, 0))
Allsim<- list()
for(i in 1:8){
  Betas<- B[[i]]
  Model<- i-1
  set.seed(0); Sim<- DetectOutbreaks::simulate(B = Betas, Model = Model, time = 60, adj.matrix = sim_adjmat)
  Allsim[[i]]<- Sim
}

#Simulated data (Figure 2 in paper)
DetectOutbreaks:::publicationfig2(Allsim)

#Fit models (If you have a GPU, use "GPU = TRUE" instead for fast computations)
Simulationfits<- list()
for(i in 1:8){
  currentsim<- Allsim[[i]]
  Model<- i-1
  HMCFit<- DetectOutbreaks::infer(y = currentsim[[1]], e_it = currentsim[[2]], Model = Model, adjmat = sim_adjmat, Stan = T, iter = 5000, GPU = F, verbose = T)
  Simulationfits[[i]]<- HMCFit
}

# Posterior fits and credible intervals (Figure 3 in paper)
DetectOutbreaks:::publicationfig3(all.infobjects=Simulationfits, all.simobjects=Allsim, adjmat=sim_adjmat)

# Log-marginal likelihoods for simulation study (Table 3 in paper)
Excessfits<- list()
for(i in 1:8){
  currentsim<- Allsim[[i]]
  index<- i*8-8
  for(j in 1:8){
    Model<- j-1
    HMCFit<- DetectOutbreaks::infer(y = currentsim[[1]], e_it = currentsim[[2]], Model = Model, adjmat = sim_adjmat, Stan = T, iter = 5000, GPU = F, verbose = T)
    Excessfits[[index+j]]<- HMCFit
  }
}
for(i in 1:8){
  indices<- i*8
  indices<- (indices-8+1):indices
  currentfits<- Excessfits[indices]
  currentsim<- Allsim[[i]]
  currentsim<- currentsim[1:2]
  DetectOutbreaks:::ModComp(allmods=currentfits, alldata=currentsim, adjmat=sim_adjmat)
}

# Heat-maps for simulation study (Figure 4 in paper)
DetectOutbreaks:::publicationfig4(all.infobjects=Simulationfits, all.simobjects=Allsim, adjmat=sim_adjmat)

# ROC-curves for simulation study (Figure 5 in paper)
DetectOutbreaks:::publicationfig5(all.infobjects=Simulationfits, all.simobjects=Allsim, adjmat=sim_adjmat)

# Section 4 - Application to invasive meningococcal disease
# Monthly reported cases (Figure 6 in paper)
data(ApplicationCounts, package = "DetectOutbreaks")
data(ApplicationPopulation, package = "DetectOutbreaks")
data(ApplicationAdjMat, package = "DetectOutbreaks")
y<- ApplicationCounts
e_it<- ApplicationPopulation
Adjmat<- ApplicationAdjMat
DetectOutbreaks::sim.plot(y, names=colnames(Adjmat))

Applicationfits<- list()
for(i in 1:8){
  Model<- i-1
  HMCFit<- DetectOutbreaks::infer(y = y, e_it = e_it, Model = Model, adjmat = Adjmat, Stan = T, iter = 5000, GPU = F, verbose = T)
  Applicationfits[[i]]<- HMCFit
}

# Heat-maps for application (Figure 7 in paper)
alldata<- list(y, e_it, colnames(Adjmat))
DetectOutbreaks:::publicationfig7(all.infobjects=Applicationfits, realdata=alldata, adjmat=Adjmat)

# Posterior median relative risks (Figure 8 in paper)
DetectOutbreaks:::publicationfig8(all.infobjects=Applicationfits)

# Posterior means for trend and seasonal components (Figure 9 in paper)
DetectOutbreaks:::publicationfig9(all.infobjects=Applicationfits)

# Table 4 - Simply extract parameter estimates from each fit-object in the list "Applicationfits". e.g., Applicationfits[[1]]$print(max_rows=400)

# Log marginal likelihoods for application (Table 5 in paper)
DetectOutbreaks:::ModComp(allmods=Applicationfits, alldata=alldata, adjmat=Adjmat)

# Posterior densities for other parameters (Figure S1 in paper)
Gammas<- c(0.1, 0.2)
params<- list(rep(0,3), c(Gammas,1.65), c(Gammas,1.65), c(Gammas,1.25,0.75), c(Gammas,0.55), c(Gammas,0.49), c(Gammas,0.35, 0.20), c(Gammas,1.65))
DetectOutbreaks:::publicationfigS1(all.infobjects=Simulationfits, all.simobjects=Allsim, all.modparams=params)

# Posterior densities for spatial components (Figure S2 in paper)
DetectOutbreaks:::publicationfigS2(all.infobjects=Simulationfits, all.simobjects=Allsim)

# Posterior predictive fits for the total case counts across all spatial locations (Figure S3 in paper)
DetectOutbreaks:::publicationfigS3(all.infobjects=Excessfits, all.simobjects=Allsim, adjmat=sim_adjmat)

# Posterior predictive fits for application (Figure S4)
DetectOutbreaks:::publicationfigS4(all.infobjects=Applicationfits,  realdata=alldata, adjmat=Adjmat)

# Correlogram for posterior probabilities of outbreak from models I-VII (Figure S5 in paper)
OutP<- list()
for(i in 1:7){
  Currentfit<- Applicationfits[[i+1]]
  Model<- i
  OutP[[i]]<- DetectOutbreaks::OutbreakProbability(y=y, e_it=e_it,inf.object=Currentfit, adjmat=Adjmat, Model=Model)
}
OutPdf<- data.frame(mod1=c(OutP[[1]]), mod2=c(OutP[[2]]), mod3=c(OutP[[3]]), mod4=c(OutP[[4]]), mod5=c(OutP[[5]]), mod6=c(OutP[[6]]), mod7=c(OutP[[7]]))
ggpairs(OutPdf, diag = list(continuous = "blankDiag"), columnLabels = c("Model I", "Model II", "Model III","Model IV","Model V","Model VI","Model VII"))
