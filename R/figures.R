
#' A function for plotting simulation figures.
#'
#' @param sim.object Resulting simulation object from either "simulate", or "simulate.old".
#'
#' @return Figures
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' sim.plot(sim.object = TRUTHS)
#'
#' set.seed(4);
#' TRUTHS2<- simulate.old(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' sim.plot(TRUTHS2)
#'
sim.plot <- function(sim.object) {
  if (is.list(sim.object)) {
    spatdata <- sim.object[[1]]
  } else {
    spatdata <- sim.object
  }
  ts_spatdata <- as.data.frame(t(spatdata))
  ts_spatdata$Time <- 1:ncol(spatdata)
  colnames(ts_spatdata) <- c(paste("u", 1:(ncol(ts_spatdata) - 1), sep = ""), "Time")
  long_data <- reshape2::melt(ts_spatdata, id.vars = "Time")

  library(ggplot2)
  a <- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time [month/year]", y = "Case counts", color = "Location") +
    guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))

  return(a)
}


#' A function for plotting mcmc diagnostic figures.
#'
#' @param inf.object Resulting fit object from either "infer", or "infer.old.
#'
#' @return Figures
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' fit<-infer(TRUTHS[[1]],TRUTHS[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=TRUE)
#' fit<- fit[[1]]  #if fit is a list containing the outbreak probabilities as its second element.
#' mcmc.plot(inf.object = fit, OutbreakProb = TRUE)
#'
#' set.seed(4);
#' TRUTHS1<- simulate.old(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' fit1<-infer.old(TRUTHS1[[1]],TRUTHS1[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=TRUE)
#' fit1<- fit1[[1]]  #fit1 is a list containing the outbreak probabilities as its second element.
#' mcmc.plot(inf.object = fit1, OutbreakProb = TRUE)
#'
mcmc.plot<- function(inf.object){

  if(is.data.frame(inf.object)){

   # par(mfrow=c(3, 3))
  #  for (i in 1:ncol(inf.object)) {
  #    hist(inf.object[-(1:2000), i], main = colnames(inf.object)[i], xlab ="", col = "white", border = "black")
   # }

    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      plot(inf.object[, i], type = "l", main = colnames(inf.object)[i], xlab ="MCMC iterations", ylab = "", col = "purple")
      grid()
    }
  }else{
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B", "state1_stationary_dist", "log_lik"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "state1_stationary_dist", "log_lik"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    colnames(inf.object) <- gsub("^1.", "", colnames(inf.object))

    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      hist(inf.object[, i], main = colnames(inf.object)[i], xlab ="", col = "white", border = "black")
    }

    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      plot(inf.object[, i], type = "l", main = colnames(inf.object)[i], xlab ="HMC iterations", ylab = "", col = "red")
      grid()
    }
  }
}

#' A function to plot inferred trend and seasonal components from Stan/MCMC output.
#'
#' @param inf.object Resulting object from either "infer" or "infer.old" functions.
#'
#' @return Figures
#' @export
#'
#' @examples set.seed(4);
#' sim_adjmat<- matrix(0, nrow = 9, ncol = 9)
#' uppertriang<- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
#' gdata::upperTriangle(sim_adjmat, byrow=TRUE)<- uppertriang
#' gdata::lowerTriangle(sim_adjmat, byrow=FALSE)<- uppertriang
#' TRUTHS<- simulate(Model = 0, time = 48, adj.matrix = sim_adjmat)
#' fit<-infer(TRUTHS[[1]],TRUTHS[[2]],Model=0,adjmat=sim_adjmat,Stan=TRUE,ModEvid=TRUE,OutbreakProb=TRUE)
#' fit<- fit[[1]]  #fit is a list containing the outbreak probabilities as its second element.
#' inf.plot(inf.object = fit)
#'
inf.plot<- function(inf.object){

  if(is.data.frame(inf.object)){
  rPosterior<- inf.object[, startsWith(colnames(inf.object), "r")]
  sPosterior<- inf.object[, startsWith(colnames(inf.object), "s")]
  inf.r<- colMeans(rPosterior)
  inf.s<- colMeans(sPosterior)
  uCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,2]
  lCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,1]
  uCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,2]
  lCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,1]
  }else{
  inf.r<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
  uCI.r<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "r")[,1,]))[,2]
  lCI.r<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "r")[,1,]))[,1]
  inf.s<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
  uCI.s<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "s")[,1,]))[,2]
  lCI.s<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "s")[,1,]))[,1]
  inf.u<- colMeans(as.data.frame(inf.object$draws(variables = "uconstrained")[,1,]))
  uCI.u<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "uconstrained")[,1,]))[,2]
  lCI.u<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "uconstrained")[,1,]))[,1]
  }

  par(mfrow=c(1,2))
  plot(0, type = "n", xlim = c(1,length(inf.r)), ylim = c(min(lCI.r, inf.r), max(uCI.r, inf.r)), ylab = "Trend component", xlab = "Time")
  polygon(c(1:length(inf.r), rev(1:length(inf.r))), c(lCI.r, rev(uCI.r)),
          col = "pink", border = NA)
  lines(1:length(inf.r), inf.r, col="red")
  grid()

    plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(min(lCI.s, inf.s), max(uCI.s, inf.s)), ylab = "Seasonal component", xlab = "Season")
    polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
            col = "pink", border = NA)
    lines(1:length(inf.s), inf.s, col="red")
    grid()
    add_legend("topright", legend="Posterior means", lty=1, col="red",
               horiz=TRUE, bty='n', cex=1.1)
}

publicationfigs0<- function(all.simobjects){
  plotlists<- list()

  for(i in 1:8){
    Model<- i
    sim.object<- all.simobjects[[i]]

    spatdata<- sim.object[[1]]

    ts_spatdata <- as.data.frame(t(spatdata))
    ts_spatdata$Time <- 1:ncol(spatdata)
    naming<- c(paste("u", 1:(ncol(ts_spatdata)-1), sep=""), "Time")
    colnames(ts_spatdata)<- naming

    Colors <- rep(c("blue", "red"), length.out = nrow(spatdata))
    Linetypes <- rep(c("dotted", "dashed", "dotdash", "longdash", "twodash"), length.out = nrow(spatdata))

    long_data <- reshape2::melt(ts_spatdata, id.vars = "Time")
    library(ggplot2)
    if(Model==8){
    rfigs<- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable, linetype = variable)) +
           geom_line() +
           scale_color_manual(values = Colors) +
           scale_linetype_manual(values = Linetypes) +
           ylim(0, 80) +
           labs(x = "Time [month]", y = "Case counts", color = "Location") +
           guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
      theme(axis.title.y = element_text(size=18),
            axis.title.x = element_text(size=18),
            axis.text.x = element_text(size=16),
            axis.text.y = element_text(size=16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16))
    }else{
      rfigs<- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable, linetype = variable)) +
        geom_line() +
        scale_color_manual(values = Colors) +
        scale_linetype_manual(values = Linetypes) +
        ylim(0, 80) +
        labs(x = "Time [month]", y = "Case counts", color = "Location") +
        guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
        theme(axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
          legend.position = "none")
    }
    plotlists[[Model]]<- rfigs
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:4], ncol = 4, labels = c("A", "B", "C", "D"), label_size = 17)
  row_2<- cowplot::plot_grid(plotlist = plotlists[5:8], ncol = 4, labels = c("E", "F", "G", "H"), label_size = 17, rel_widths = c(0.75, 0.75, 0.75, 1))
  cowplot::plot_grid(row_1, row_2, nrow = 2)
}


publicationfigs1<- function(all.infobjects, all.simobjects, adjmat){
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])

  pdf("AllfitsB.pdf", paper="a4", width=12,height=12, pointsize=12)
  par(mfrow=c(4,3))

  for(i in c(2,4,6,8)){
    sim.object<- all.simobjects[[i]]
    inf.object<- all.infobjects[[i]]
    Model<- i-1

    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    sim.r<- sim.object[[3]]
    sim.s<- sim.object[[4]]
    sim.u<- sim.object[[5]]
    X_it<- sim.object[[6]]

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

    sum_Y<- matrix(NA, nrow = length(thinning), ncol = time)
    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    if(Model == 0){
      for(index in 1:length(thinning)){
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Xit<- sum_Xit + Ex_Xit
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Y[index, ]<- colSums(pred_Y)
      }
    }
    inf.r<- colMeans(fullr.draws)
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(fullr.draws))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(fullr.draws))[,1]
    inf.s<- colMeans(fulls.draws)
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(fulls.draws))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(fulls.draws))[,1]
    inf.u<- colMeans(fullu.draws)
    uCI.u<- posterior_interval_custom(as.matrix.data.frame(fullu.draws))[,2]
    lCI.u<- posterior_interval_custom(as.matrix.data.frame(fullu.draws))[,1]
    inf.Y<- colMeans(sum_Y)
    uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,2]
    lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,1]

    plot(0, type = "n", xlim = c(1,length(sim.r)), ylim = c(-14.5, -12.8), ylab = "trend component", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7)
    polygon(c(1:length(sim.r), rev(1:length(sim.r))), c(lCI.r, rev(uCI.r)),
            col = "pink", border = NA)
    lines(1:length(inf.r), inf.r, col="red", lty=1)
    points(1:length(sim.r), sim.r, pch = 19)
    grid()
    legendary::labelFig(LETTERS[Model+1], adj = c(-0.15, 0.10), font=2, cex=1.8)

    if(length(inf.s) < length(sim.s)){
      plot(0, type = "n", xlim = c(1,length(sim.s)), ylim = c(-2.0, 1.5), ylab = "seasonal component", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7)
      polygon(c(1:length(sim.s), rev(1:length(sim.s))), c(rep(lCI.s,round(length(sim.s)/length(inf.s))), rev(rep(uCI.s, round(length(sim.s)/length(inf.s))))),
              col = "pink", border = NA)
      lines(1:length(sim.s), rep(inf.s, round(length(sim.s)/length(inf.s))), col="red", lty=1)
    }else if(length(inf.s) > length(sim.s)){
      plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(-2.0, 1.5), ylab = "seasonal component", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7)
      polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
              col = "pink", border = NA)
      lines(1:length(inf.s), inf.s, col="red", lty=1)
    }else{
      plot(0, type = "n", xlim = c(1,length(sim.s)), ylim = c(-2.0, 1.5), ylab = "seasonal component", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7)
      polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
              col = "pink", border = NA)
      lines(1:length(sim.s), inf.s, col="red", lty=1)
    }
    points(1:length(sim.s), sim.s, pch = 19)
    grid()

    plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(0, 300), ylab = "overall case counts", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7)
    polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
    points(1:ncol(y), colSums(y), pch = 19)
    grid()
  }
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}


publicationfigs2<- function(all.infobjects, all.simobjects){
  plotlists<- list()

  for(i in 1:8){
    Model<- i-1
    inf.object<- all.infobjects[[i]]
    sim.object<- all.simobjects[[i]]

    if(!is.data.frame(inf.object)){
      fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
    }else{
      fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
    }

    thinning<- numeric(floor(nrow(fullu.draws)/10))
    thinning[1]<- 10
    for(i in 2:length(thinning)){
      thinning[i]<- thinning[i-1] + 10
    }
    u.draws<- fullu.draws[thinning, ]
    sim.u<- sim.object[[5]]

  # Violin plot for spatial components
  spatcomp<- data.frame(value = c(u.draws[, 1], u.draws[, 2], u.draws[, 3], u.draws[, 4],
                                  u.draws[, 5], u.draws[, 6], u.draws[, 7], u.draws[, 8],
                                  u.draws[, 9]), group = factor(rep(c("u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8", "u9"), each = nrow(u.draws))))
  spatcomp$group <- factor(spatcomp$group, levels = unique(spatcomp$group))
  library(RColorBrewer)
  mycolors <- c(rep(c("blue", "red"), 4), "blue")
  rfigs<- ggplot(spatcomp, aes(x = group, y = value, fill = group)) +
          geom_violin(trim = FALSE, alpha = 0.7) +
          geom_point(x = 1, y = sim.u[1], size = 2, shape = 19) +
          geom_point(x = 2, y = sim.u[2], size = 2, shape = 19) +
          geom_point(x = 3, y = sim.u[3], size = 2, shape = 19) +
          geom_point(x = 4, y = sim.u[4], size = 2, shape = 19) +
          geom_point(x = 5, y = sim.u[5], size = 2, shape = 19) +
          geom_point(x = 6, y = sim.u[6], size = 2, shape = 19) +
          geom_point(x = 7, y = sim.u[7], size = 2, shape = 19) +
          geom_point(x = 8, y = sim.u[8], size = 2, shape = 19) +
          geom_point(x = 9, y = sim.u[9], size = 2, shape = 19) +
          ylim(-0.90, 0.70) +
          labs(title = "", x = "Location", y = "Value", fill = "") +
          theme_minimal() +
          scale_fill_manual(values = mycolors) +
          theme(axis.title.y = element_text(size=18),
                axis.title.x = element_text(size=18),
                axis.text.x = element_text(size=16),
                axis.text.y = element_text(size=16),
                legend.title = element_text(size = 18),
                legend.text = element_text(size = 16),legend.position = "none")
  plotlists[[Model+1]]<- rfigs
  }
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 4, labels = c("A", "B", "C", "D", "E", "F", "G", "H"), label_size = 17))
  add_legend(0.85, 1.20, legend="Truth",
             pch=19, col="black",
             horiz=TRUE, bty='n', cex=1.5)
  #add_legend("topright", legend="Truth",
  #           pch=19, col="black",
  #           horiz=TRUE, bty='n', cex=1.1)
}

publicationfigs3<- function(all.infobjects, all.simobjects, all.modparams){
  library(ggplot2)
  plotlists<- list()

  for(i in 1:7){
    Model<- i
    inf.object<- all.infobjects[[i+1]]
    sim.object<- all.simobjects[[i+1]]
    modparams<- all.modparams[[i+1]]

    if(!is.data.frame(inf.object)){
      fullG12.draws<- stack(as.data.frame(inf.object$draws(variables = "G12")[,1,]))[,1]
      fullG21.draws<- stack(as.data.frame(inf.object$draws(variables = "G21")[,1,]))[,1]
      fulld.draws<- stack(as.data.frame(inf.object$draws(variables = "state1_stationary_dist")[,1,]))[,1]
    }else{
      fullG12.draws<- as.numeric(inf.object[-(1:burn.in), 1])
      fullG21.draws<- as.numeric(inf.object[-(1:burn.in), 2])
      fulld.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+2+1]
    }

    thinning<- numeric(floor(length(fullG12.draws)/10))
    thinning[1]<- 10
    for(i in 2:length(thinning)){
      thinning[i]<- thinning[i-1] + 10
    }

    G12.draws<- fullG12.draws[thinning]
    G21.draws<- fullG21.draws[thinning]
    d.draws<- fulld.draws[thinning]

    if(Model %in% c(1,2,4,5,7)){
      if(!is.data.frame(inf.object)){
        B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
        B.draws<- B.draws[thinning]
      }else{
        B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
        B.draws<- B.draws[thinning]
      }
    }else if(Model %in% c(3,6)){
      if(!is.data.frame(inf.object)){
        B.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
        B.draws<- B.draws[thinning, ]
      }else{
        B.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+(1:2)]
        B.draws<- B.draws[thinning, ]
      }
    }

    # Violin plot
    if(Model %in% c(1,2,4,5,7)){
      vioB<- data.frame(value = c(G12.draws, G21.draws, d.draws, B.draws), group = factor(rep(c("G12", "G21", "d1", "B1"), each = length(B.draws))))
      vioB$group <- factor(vioB$group, levels = unique(vioB$group))
      dTr<- state_dist_cpp(modparams[1], modparams[2])[2]
      vioBTr<- modparams[3]
      library(RColorBrewer)
      rfigs<- ggplot(vioB, aes(x = group, y = value, fill = group)) +
              geom_violin(trim = FALSE, alpha = 0.7) +
              geom_point(x = 1, y = modparams[1], size = 2, shape = 19) +
              geom_point(x = 2, y = modparams[2], size = 2, shape = 19) +
              geom_point(x = 3, y = dTr, size = 2, shape = 19) +
              geom_point(x = 4, y = vioBTr, size = 2, shape = 19) +
              ylim(0.0, 2.0) +
              labs(title = "", x = "Parameter", y = "Value", fill = "True value") +
              theme_minimal() +
              scale_fill_brewer(palette = "Set2") +
              scale_x_discrete(labels = c("G12" = expression(gamma[0][1]),
                                          "G21" = expression(gamma[1][0]),
                                          "d1" = expression(delta[1]),
                                          "B1" = expression(beta[1]))) +
              scale_fill_manual(values = c(rep("purple", 4)),
                                labels = c(expression(gamma[0][1]),
                                           expression(gamma[1][0]),
                                           expression(delta[1]),
                                           expression(beta[1])),
                                name = "True value")+
              theme(axis.title.y = element_text(size=18),
                    axis.title.x = element_text(size=18),
                    axis.text.x = element_text(size=16),
                    axis.text.y = element_text(size=16),
                    legend.title = element_text(size = 18),
                    legend.text = element_text(size = 16), legend.position = "none")
    }else{
      vioB<- data.frame(value = c(G12.draws, G21.draws, d.draws, B.draws[,1], B.draws[,2]), group = factor(rep(c("G12", "G21", "d1", "B1", "B2"), each = length(G12.draws))))
      vioB$group <- factor(vioB$group, levels = unique(vioB$group))
      dTr<- state_dist_cpp(modparams[1], modparams[2])[2]
      vioBTr<- c(modparams[3], modparams[4])
      library(RColorBrewer)
      rfigs<- ggplot(vioB, aes(x = group, y = value, fill = group)) +
              geom_violin(trim = FALSE, alpha = 0.7) +
              geom_point(x = 1, y = modparams[1], size = 2, shape = 19) +
              geom_point(x = 2, y = modparams[2], size = 2, shape = 19) +
              geom_point(x = 3, y = dTr, size = 2, shape = 19) +
              geom_point(x = 4, y = vioBTr[1], size = 2, shape = 19) +
              geom_point(x = 5, y = vioBTr[2], size = 2, shape = 19) +
              ylim(0.0, 2.0) +
              labs(title = "", x = "Parameter", y = "Value", fill = "True value") +
              theme_minimal() +
              scale_fill_brewer(palette = "Set2") +
              scale_x_discrete(labels = c("G12" = expression(gamma[0][1]),
                                          "G21" = expression(gamma[1][0]),
                                          "d1" = expression(delta[1]),
                                          "B1" = expression(beta[1]),
                                          "B2" = expression(beta[2]))) +
              scale_fill_manual(values = c(rep("purple", 5)),
                                labels = c(expression(gamma[0][1]),
                                           expression(gamma[1][0]),
                                           expression(delta[1]),
                                           expression(beta[1]),
                                           expression(beta[2])),
                                name = "True value") +
              theme(axis.title.y = element_text(size=18),
                    axis.title.x = element_text(size=18),
                    axis.text.x = element_text(size=16),
                    axis.text.y = element_text(size=16),
                    legend.title = element_text(size = 18),
                    legend.text = element_text(size = 16), legend.position = "none")
    }
    plotlists[[Model]]<- rfigs
  }
  plot.new()
  print(cowplot::plot_grid(plotlist = plotlists[1:7], ncol = 4, labels = c("A", "B", "C", "D","E","F","G"), label_size = 17))
  #row_1<- cowplot::plot_grid(plotlist = plotlists[1:4], ncol = 4, labels = c("A", "B", "C", "D"))#, rel_widths = c(0.95, 0.95, 1, 0.95))
  #row_2<- cowplot::plot_grid(plotlist = plotlists[5:7], ncol = 3, labels = c("E", "F", "G"))#, rel_widths = c(0.95, 1, 0.95))
  #print(cowplot::plot_grid(row_1, row_2, nrow = 2))
  add_legend(0.73, -0.4, legend="Truth",
             pch=19, col="black",
             horiz=TRUE, bty='n', cex=1.8)
  #legend(0.75, -0.4, legend = c("Truth"), col = "black", pch = 19, cex=1.2)
}


publicationfigs4<- function(all.infobjects, all.simobjects, adjmat){
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])
  #pdf("posterheatmap.pdf", paper="special", width=21,height=12, pointsize=14)
  plot.new()
  par(mfrow=c(2,4))
  X_it<- all.simobjects[[2]][[6]]
  smallxit<- X_it[c(1,3,5,7,9), ]
  bigxit<- X_it[c(2,4,6,8), ]
  bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]
  image(x=1:time, y=1:ndept, t(bigsmallxit), main ="", axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80)

  abline(h=4.5, col="black", lty=2)
  #custom Y-axis
  axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
  axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
  #custom X-axis
  axis(1, cex.axis = 1.8)
  legendary::labelFig("A", adj = c(-0.15, 0.10), font=2, cex=1.8)

  for(i in 1:7){
    sim.object<- all.simobjects[[i+1]]
    inf.object<- all.infobjects[[i+1]]
    Model<- i

  y<- sim.object[[1]]
  e_it<- sim.object[[2]]
  sim.r<- sim.object[[3]]
  sim.s<- sim.object[[4]]
  sim.u<- sim.object[[5]]
  X_it<- sim.object[[6]]

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

  sum_Xit<- matrix(0, nrow = ndept, ncol = time)

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
      Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B.draws[index], Model = Model, adjmat = adjmat, Cyclic = T)
      sum_Xit<- sum_Xit + Ex_Xit
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
    }
  }

  if(Model != 0){
    mean_Xit<- sum_Xit/length(thinning)
    smallmeanxit<- mean_Xit[c(1,3,5,7,9), ]
    bigmeanxit<- mean_Xit[c(2,4,6,8), ]
    bigsmallmeanxit<- mean_Xit[c(2,4,6,8,1,3,5,7,9), ]
    image(x=1:time, y=1:ndept, t(bigsmallmeanxit), main = "", axes=F, ylab = "spatial location", xlab = "Time [month]", cex.lab=1.80)
    abline(h=4.5, col="black", lty=2)
    #custom Y-axis
    axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    #custom X-axis
    axis(1, cex.axis = 1.8)
    legendary::labelFig(LETTERS[Model+1], adj = c(-0.15, 0.10), font=2, cex=1.8)
  }
  }
  #dev.off()
}

publicationfigs5<- function(all.infobjects, all.simobjects, adjmat){
  plotlists<- list()
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])
  par(mfrow=c(2,4))
  X_it<- all.simobjects[[2]][[6]]
  smallxit<- X_it[c(1,3,5,7,9), ]
  bigxit<- X_it[c(2,4,6,8), ]
  bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]

  for(i in 1:7){
    sim.object<- all.simobjects[[i+1]]
    inf.object<- all.infobjects[[i+1]]
    Model<- i

    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    sim.r<- sim.object[[3]]
    sim.s<- sim.object[[4]]
    sim.u<- sim.object[[5]]
    X_it<- sim.object[[6]]

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

    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

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
        Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B.draws[index], Model = Model, adjmat = adjmat, Cyclic = T)
        sum_Xit<- sum_Xit + Ex_Xit
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
      }
    }

    if(Model != 0){
      #par(mfrow=c(2,4))
      mean_Xit<- sum_Xit/length(thinning)
      smallmeanxit<- mean_Xit[c(1,3,5,7,9), ]
      bigmeanxit<- mean_Xit[c(2,4,6,8), ]
      bigsmallmeanxit<- mean_Xit[c(2,4,6,8,1,3,5,7,9), ]
      # Load required libraries
      library(pROC)
      library(ggplot2)
      library(dplyr)

      # Create ROC data
      smallroc_data <- pROC::roc(as.vector(smallxit), as.vector(smallmeanxit))
      bigroc_data <- pROC::roc(as.vector(bigxit), as.vector(bigmeanxit))

      # Extract ROC data into data frames
      smallroc_df <- data.frame(
        FPR = 1 - smallroc_data$specificities,
        TPR = smallroc_data$sensitivities,
        Group = "Small Cities"
      )
      bigroc_df <- data.frame(
        FPR = 1 - bigroc_data$specificities,
        TPR = bigroc_data$sensitivities,
        Group = "Large Cities"
      )

      # Combine both data frames
      roc_df <- bind_rows(smallroc_df, bigroc_df)

      # Calculate AUC values
      auc_small <- round(auc(smallroc_data), 2)
      auc_big <- round(auc(bigroc_data), 2)

      # Plot ROC curves with ggplot2
      rfig<- ggplot(roc_df, aes(x = FPR, y = TPR, color = Group)) +
        geom_line(size = 1.1) +
        scale_color_manual(values = c("Small Cities" = "blue", "Large Cities" = "red")) +
        labs(
          x = "1 - Specificity",
          y = "Sensitivity",
          title = "",
          color = "City Size"
        ) +
        annotate("text", x = 0.60, y = 0.10,
                 label = paste("AUC small cities =", auc_small),
                 color = "blue") +
        annotate("text", x = 0.60, y = 0.05,
                 label = paste("AUC large cities =", auc_big),
                 color = "red") +
        theme_minimal() +
        theme(axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 16), legend.position = "bottomright") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")

      plotlists[[Model]]<- rfig
    }
  }
  cowplot::plot_grid(plotlist = plotlists, ncol = 4, labels = c("A", "B", "C", "D", "E", "F", "G"), label_size = 17)
}

publicationfigs6<- function(all.infobjects,  realdata, adjmat){
  time<- ncol(realdata[[1]])
  ndept<- nrow(realdata[[1]])

  pdf("Correctedpostpred.pdf", paper="special", width=18,height=9, pointsize=12)
  par(mfrow=c(2,4))

  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i-1

    y<- realdata[[1]]
    e_it<- realdata[[2]]

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

    sum_Y<- matrix(NA, nrow = length(thinning), ncol = time)
    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    if(Model == 0){
      for(index in 1:length(thinning)){
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Xit<- sum_Xit + Ex_Xit
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Y[index, ]<- colSums(pred_Y)
      }
    }
    inf.r<- colMeans(fullr.draws)
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(fullr.draws))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(fullr.draws))[,1]
    inf.s<- colMeans(fulls.draws)
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(fulls.draws))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(fulls.draws))[,1]
    inf.u<- colMeans(fullu.draws)
    uCI.u<- posterior_interval_custom(as.matrix.data.frame(fullu.draws))[,2]
    lCI.u<- posterior_interval_custom(as.matrix.data.frame(fullu.draws))[,1]
    inf.Y<- colMeans(sum_Y)
    uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,2]
    lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,1]

    plot(0, type = "n", xaxt = "n", xlim = c(1,ncol(y)), ylim = c(min(lCI.Y, y, na.rm = T), max(uCI.Y, y, na.rm = T)), ylab = "overall case counts", xlab = "Time [month/year]", cex.axis = 1.8, cex.lab=1.8)
    polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
    points(1:ncol(y), colSums(y, na.rm = T), pch = 19)
    years<- 2013:2019
    axis(1, at = seq(1, 84, by = 12), labels = years, cex.axis = 1.8)
    grid()
    legendary::labelFig(LETTERS[Model+1], adj = c(-0.15, 0.10), font=2, cex=1.8)
  }
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}

publicationfigs7<- function(all.infobjects, realdata, adjmat){
  y<- realdata[[1]]
  e_it<- realdata[[2]]
  countries<- realdata[[3]]
  time<- ncol(y)
  ndept<- nrow(y)
  listofOutProbs<- list()

  #pdf("Allconnectedheatmap.pdf", paper="special", width=21,height=12, pointsize=14)
  #par(mfrow=c(2,4))

    for(i in 1:7){
    inf.object<- all.infobjects[[i+1]]
    Model<- i

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

    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

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
        Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B.draws[index], Model = Model, adjmat = adjmat, Cyclic = T)
        sum_Xit<- sum_Xit + Ex_Xit
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
      }
    }


    #if(Model != 0){
    #  mean_Xit<- sum_Xit/length(thinning)
    #  mean_Xit<- mean_Xit[ndept:1, ]
    #  par(mar = c(4, 7.5, 4, 1))
    #  image(x=1:time, y=1:ndept, t(mean_Xit), main = "", axes=F, ylab = "", xlab = "Time [month/year]", cex.lab=1.8)
      #custom Y-axis
    #  axis(2, at=seq(1, length(countries), length.out=length(countries)), labels=rev(countries), lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.40)
      #custom X-axis
    #  years<- 2013:2019
    #  axis(1, at = seq(1, 84, by = 12), labels = years, cex.axis = 1.8)
    #  legendary::labelFig(LETTERS[Model], adj = c(-0.15, 0.05), font=2, cex=2.0)
    #}
    listofOutProbs[[Model]]<- sum_Xit/length(thinning)
    }
  #dev.off()
  return(listofOutProbs)
}


publicationfigs8<- function(all.infobjects){
  plotlists<- list()

  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i

    if(!is.data.frame(inf.object)){
      fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
    }else{
      fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
    }

#get required regions
  poly <- cshapes::cshp(date=as.Date("2019-12-31"), useGW=TRUE)
  #newpoly<- poly[poly$country_name %in% c("Austria","Belgium","Cyprus","Czech Republic",
  #                                      "Denmark","Estonia","Finland","France","German Federal Republic",
  #                                      "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
  #                                      "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
  #                                      "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom"), ]


  newpoly<- poly[poly$country_name %in% c("Albania", "Austria", "Belarus (Byelorussia)", "Belgium", "Bosnia-Herzegovina",
                                          "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland",
                                          "France", "German Federal Republic", "Greece", "Hungary", "Iceland",
                                          "Ireland", "Italy/Sardinia", "Kosovo", "Latvia", "Lithuania", "Luxembourg", "Malta",
                                          "Macedonia (FYROM/North Macedonia)", "Moldova", "Montenegro", "Netherlands", "Norway",
                                          "Poland", "Portugal", "Rumania", "Russia (Soviet Union)", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
                                          "Switzerland", "Ukraine", "United Kingdom"), ]

  Europe<- c("Austria","Belgium","Cyprus","Czech Republic",
             "Denmark","Estonia","Finland","France","German Federal Republic",
             "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
             "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
             "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom",
             "Albania", "Belarus (Byelorussia)", "Bosnia-Herzegovina", "Bulgaria",
             "Croatia", "Kosovo", "Russia (Soviet Union)",
             "Moldova", "Montenegro", "Macedonia (FYROM/North Macedonia)",
             "Serbia", "Switzerland", "Ukraine")

  sortedpoly <- newpoly[order(newpoly$country_name), ]
  #sortedpoly$country_name<- c("Austria","Belgium","Cyprus","Czechia",
  #                          "Denmark","Estonia","Finland","France","Germany",
  #                          "Greece","Hungary","Iceland","Ireland","Italy","Latvia",
  #                          "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
  #                          "Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","United Kingdom")

  u.draws<- colMeans(fullu.draws)
  #expUi<- c(exp(u.draws), rep(NA, 13))
  expUi<- c(u.draws, rep(NA, 13))
  uidf<- data.frame(expUi = expUi, countryname = Europe)
  sorteduidf<- uidf[order(uidf$countryname), ]
  rrDF<- data.frame(rr=sorteduidf$expUi, gwcode=sortedpoly$gwcode, index=c(1:nrow(uidf)), name=sorteduidf$countryname)

  library(sf)
  sortedpoly<- st_make_valid(newpoly)
  bbox <- st_bbox(c(xmin = -29, ymin = 30, xmax = 40.17875, ymax = 71.15471), crs = st_crs(sortedpoly))
  sortedpoly<- st_crop(sortedpoly, bbox)
  shapefile <- ggplot2::fortify(sortedpoly, region = 'gwcode')

  shp_rrDF<- sp::merge(shapefile, rrDF,
                     by.x="gwcode",
                     by.y="gwcode",
                     all.x=F)

  library(ggplot2)
  library(viridis)
  if(Model == 8){
   rfig<- (ggplot(data = shp_rrDF) +
              geom_sf(aes(fill = rr)) +
              #geom_sf_text(aes(label = country_name), size = 3) +
              #scale_fill_viridis(option = "turbo", direction = 1, alpha=1, begin=0.6, end=1, na.value = "lightgrey") +  # Reverse the color scale
              scale_fill_gradient2(low = "lightblue", mid = "blue", high = "red", midpoint = 0, na.value = "lightgrey", labels = ~ format(round(exp(.),1), nsmall=1)) +
              coord_sf() + theme_void() +
              ggtitle("") + theme(plot.title = element_text(hjust = 0.55)) +
              labs(fill = expression(Exp(u[i]))) +
             theme(legend.title = element_text(size = 18),
                   legend.text = element_text(size = 16)))
  }else{
  rfig<- (ggplot(data = shp_rrDF) +
    geom_sf(aes(fill = rr)) +
    #geom_sf_text(aes(label = country_name), size = 3) +
    #scale_fill_viridis(option = "turbo", direction = 1, alpha=1, begin=0.6, end=1, na.value = "lightgrey") +  # Reverse the color scale
    scale_fill_gradient2(low = "lightblue", mid = "blue", high = "red", midpoint = 0, na.value = "lightgrey") +
    coord_sf() + theme_void() +
    ggtitle("") + theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill = "Median relative risk") +
    theme(legend.position = "none")) #+
    #theme_minimal()
  }
  plotlists[[Model]]<- rfig
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:4], ncol = 4, labels = c("A", "B", "C", "D"), label_size = 17)
  row_2<- cowplot::plot_grid(plotlist = plotlists[5:8], ncol = 4, labels = c("E", "F", "G", "H"), label_size = 17, rel_widths = c(1.16, 1.16, 1.16, 1.60))
  print(cowplot::plot_grid(row_1, row_2, nrow = 2))
}

publicationfigs9<- function(all.infobjects){
  rmat<- matrix(NA, nrow = 8, ncol = 84)
  smat<- matrix(NA, nrow = 8, ncol = 12)
  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i

    if(!is.data.frame(inf.object)){
      r.draws<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
      s.draws<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
    }else{
      r.draws<- colMeans(inf.object[-(1:burn.in), 5+(1:time)])
      s.draws<- colMeans(inf.object[-(1:burn.in), 5+time+(1:12)])
    }
    rmat[i, ]<- r.draws
    smat[i, ]<- s.draws
  }
    r_names<- c("0", "I", "II", "III", "IV", "V", "VI", "VII")
    ts_rdata <- as.data.frame(t(rmat))
    ts_sdata <- as.data.frame(t(smat))
    ts_rdata$Time <- seq.Date(from = as.Date("2013-01-01"), to = as.Date("2019-12-01"), by = "month")
    ts_sdata$Time <- seq.Date(from = as.Date("2019-01-01"), to = as.Date("2019-12-01"), by = "month")

    colnames(ts_rdata) <- c(r_names, "Time")
    colnames(ts_sdata) <- c(r_names, "Time")

    long_data <- reshape2::melt(ts_rdata, id.vars = "Time")
    library(ggplot2)
    a<- (ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
           geom_line() +
           labs(x = "Time [month/year]", y = "Trend component", color = "Model") +
           scale_x_date(date_labels = "%b %Y", date_breaks = "1 years")+
           theme(axis.title.y = element_text(size=18),
                 axis.title.x = element_text(size=18),
                 axis.text.x = element_text(size=16),
                 axis.text.y = element_text(size=16),
                 legend.title = element_text(size = 18),
                 legend.text = element_text(size = 16), legend.position = "none"))

    long_data2 <- reshape2::melt(ts_sdata, id.vars = "Time")
    library(ggplot2)
    b<- (ggplot2::ggplot(data = long_data2, mapping = aes(x = Time, y = value, color = variable)) +
           geom_line() +
           labs(x = "Time [month]", y = "Seasonal component", color = "Model") +
           scale_x_date(date_labels = "%b", date_breaks = "1 months") +
           theme(axis.title.y = element_text(size=18),
                 axis.title.x = element_text(size=18),
                 axis.text.x = element_text(size=16),
                 axis.text.y = element_text(size=16),
                 legend.title = element_text(size = 18),
                 legend.text = element_text(size = 16)))
    plotlists<- list(a, b)
    print(cowplot::plot_grid(plotlist = plotlists, ncol = 2, labels = c("A", "B"), label_size = 17, rel_widths = c(1.08, 1.15)))
}

publicationfigs10<- function(all.infobjects, all.simobjects, adjmat){
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])

  pdf("supplmA.pdf", paper="special", width=24,height=24, pointsize=12)
  par(mfrow=c(8,8))

  for(j in 1:8){
    sim.object<- all.simobjects[[j]]
    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    ind<- j*8-8

  for(k in 1:8){
    inf.object<- all.infobjects[[ind+k]]
    Model<- k-1

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

    thinning<- numeric(floor(nrow(fullr.draws)/100))
    thinning[1]<- 100
    for(i in 2:length(thinning)){
      thinning[i]<- thinning[i-1] + 100
    }

    G12.draws<- fullG12.draws[thinning]
    G21.draws<- fullG21.draws[thinning]
    r.draws<- fullr.draws[thinning, ]
    s.draws<- fulls.draws[thinning, ]
    u.draws<- fullu.draws[thinning, ]
    d.draws<- fulld.draws[thinning]

    sum_Y<- matrix(NA, nrow = length(thinning), ncol = time)
    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    if(Model == 0){
      for(index in 1:length(thinning)){
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Xit<- sum_Xit + Ex_Xit
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Y[index, ]<- colSums(pred_Y)
      }
    }
    inf.Y<- colMeans(sum_Y)
    uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,2]
    lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,1]

    plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(0, 300), ylab = "overall case counts", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7)
    polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
    points(1:ncol(y), colSums(y), pch = 19)
    grid()
  }
    legendary::labelFig(LETTERS[j], adj = c(-9.65, 0.10), font=2, cex=2.0)
  }
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}

crudevsfittedRS<- function(alldata, all.infobjects){
  y<- alldata[[1]]
  e_it<- alldata[[2]]
  crudeRes<- crudeEst(y = y, e_it = e_it)
  crudeR<- crudeRes[[1]]
  crudeS<- crudeRes[[2]]
  rmat<- matrix(NA, nrow = 9, ncol = ncol(y))
  smat<- matrix(NA, nrow = 9, ncol = 12)
  rmat[1, ]<- crudeR
  smat[1, ]<- crudeS[1:12]
  for(i in 2:9){
    inf.object<- all.infobjects[[i-1]]
    Model<- i

    if(!is.data.frame(inf.object)){
      r.draws1<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
      r.draws2<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,2,]))
      r.draws3<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,3,]))
      r.draws4<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,4,]))
      s.draws1<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
      s.draws2<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,2,]))
      s.draws3<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,3,]))
      s.draws4<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,4,]))
      r.draws<- (r.draws1 + r.draws2 + r.draws3 + r.draws4)/4
      s.draws<- (s.draws1 + s.draws2 + s.draws3 + s.draws4)/4

      #r.draws<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
      #s.draws<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
    }else{
      r.draws<- colMeans(inf.object[-(1:burn.in), 5+(1:time)])
      s.draws<- colMeans(inf.object[-(1:burn.in), 5+time+(1:12)])
    }
    rmat[i, ]<- r.draws
    smat[i, ]<- s.draws
  }
  r_names<- c("Crude", "0", "I", "II", "III", "IV", "V", "VI", "VII")
  ts_rdata <- as.data.frame(t(rmat))
  ts_sdata <- as.data.frame(t(smat))
  ts_rdata$Time <- seq.Date(from = as.Date("1999-01-01"), to = as.Date("2019-12-01"), by = "month")
  ts_sdata$Time <- seq.Date(from = as.Date("2019-01-01"), to = as.Date("2019-12-01"), by = "month")

  colnames(ts_rdata) <- c(r_names, "Time")
  colnames(ts_sdata) <- c(r_names, "Time")

  long_data <- reshape2::melt(ts_rdata, id.vars = "Time")
  library(ggplot2)
  a<- (ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
         geom_line() +
         labs(x = "Time [month/year]", y = "Trend component", color = "Model") +
         scale_x_date(date_labels = "%b %Y", date_breaks = "3 years")+
         theme(legend.position = "none"))

  long_data2 <- reshape2::melt(ts_sdata, id.vars = "Time")
  library(ggplot2)
  b<- (ggplot2::ggplot(data = long_data2, mapping = aes(x = Time, y = value, color = variable)) +
         geom_line() +
         labs(x = "Time [month]", y = "Seasonal component", color = "Model") +
         scale_x_date(date_labels = "%b", date_breaks = "1 months"))
  plotlists<- list(a, b)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 2, labels = c("A", "B"), rel_widths = c(1, 1.15)))
}
