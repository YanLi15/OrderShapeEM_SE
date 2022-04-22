## Plot the FDR and power of different methods for real data-guided simulations

n.rep    = 100
n.gene   = 10000
alphas   = seq(0, 0.1, 0.005)
sig.str  = "Moderate"
n.signal = n.gene*c(0.075, 0.1, 0.15)  # Low: 7.5%; Medium: 10%; High: 15%
info.str = c("Uninformative", "Weak", "Strong") # for single auxiliary covariate sequence

res.obj  = list()
methods <- c("BY", "BH","OrderShapeEM", "SABHA")
for (i in 1:length(info.str)){
  for (j in 1:length(n.signal)){
    data.name <- paste0("./simulation/",info.str[i], "_Sig", n.signal[j]/n.gene,"_N",n.gene,".RData")
    load(file.path(data.name))
    
    res.obj$fdr[[info.str[i]]][[paste0(n.signal[j]*100/n.gene,"%")]] <- data.frame(BY = res$BY$fdr, 
                                                                                   BH = res$BH$fdr, 
                                                                                   SABHA = res$SABHA$fdr, 
                                                                                   OrderShapeEM = res$OrderShapeEM$fdr)
    res.obj$power[[info.str[i]]][[paste0(n.signal[j]*100/n.gene,"%")]] <- data.frame(BY = res$BY$pwr, 
                                                                                     BH = res$BH$pwr, 
                                                                                     SABHA = res$SABHA$pwr, 
                                                                                     OrderShapeEM = res$OrderShapeEM$pwr)
  }
}

## plot FDR
pic <- paste0("./Simulation/plots/FDR_N", n.gene,".png")
png(pic, width=1000*3,height=1000*2,res=300)
par(mfrow=c(3,3), mgp = c(3,1,0), mai = c(0.4,0.4,0.3,0.3))
for (i in 1:length(info.str)) {
  for (j in 1:length(n.signal)) {
    plot(
      alphas[-1],
      res.obj$fdr[[i]][[j]]$BH[-1],
      xlab = "Target FDR level",
      ylab = "FDR",
      xlim = c(0, 0.1),
      ylim = c(0, 0.11),
      type = "l",
      lwd  = 1.5,
      col  = "cyan",
      main = paste0(info.str[i], " informativeness, ", n.signal[j]*100/n.gene, "% Signal")
    )
    lines(alphas[-1], res.obj$fdr[[i]][[j]]$BY[-1], lwd = 1.5, col = "tan")
    lines(alphas[-1], res.obj$fdr[[i]][[j]]$SABHA[-1], lwd = 1.5, col = "blue")
    lines(alphas[-1], res.obj$fdr[[i]][[j]]$OrderShapeEM[-1], lwd = 1.5, col = "red")
    lines(alphas, alphas, lwd = 1.5, lty = 2)
    
    if (i*j == 1)
      legend("topleft",
             methods[-6],
             col = c("tan", "cyan", "red", "blue"),
             lty = 1,
             lwd = 1.5)
  }
}
dev.off()

## plot power
pic <- paste0("./Simulation/plots/Power_N", n.gene,".png")
png(pic, width=1000*3,height=1000*2,res=300)
par(mfrow=c(3,3), mgp = c(3,1,0), mai = c(0.4,0.4,0.3,0.3))
for (i in 1:length(info.str)){
  for (j in 1:length(n.signal)){
    plot(alphas[-1], res.obj$power[[i]][[j]]$BH[-1],
         xlab = "Target FDR level",
         ylab = "Power",
         xlim = c(0,0.1),
         ylim = c(0.15,1),
         type = "l",
         lwd  = 1.5,
         col  = "cyan",
         main = paste0(info.str[i]," informativeness, ", n.signal[j]*100/n.gene, "% Signal"))
    lines(alphas[-1], res.obj$power[[i]][[j]]$BY[-1], lwd = 1.5, col = "tan")
    lines(alphas[-1], res.obj$power[[i]][[j]]$SABHA[-1], lwd = 1.5, col = "blue")
    lines(alphas[-1], res.obj$power[[i]][[j]]$OrderShapeEM[-1], lwd = 1.5, col = "red")
    
    if (i*j == 6)
      legend("bottomright",
             methods[-6],
             col = c("tan", "cyan", "red", "blue"),
             lty = 1,
             lwd = 1.5)
    
  }
}
dev.off()