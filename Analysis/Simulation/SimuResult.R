## Generate synthetic data guided by real data
## Compute FDR and Power of all methods across 100 replicates in different scenarios

library(OrderShapeEM)
library(SPARK)
library(compcodeR)
library(DESeq2)

source('./R/All_q_est_functions.R')
source('./R/SimuFunc.R')

n.rep    = 100
tau      = 0.35 
ipt      = 2   
n.gene   = 10000
sig.str  = "Moderate"     
n.signal = 1500             # 7.5% signal: 750; 10% signal: 1000; 15% signal: 1500
n.sample = 50
info.str = "Medium"  # "Uninformative", "Weak", "Medium", "Strong"
alphas   = seq(0, 0.1, 0.005)

methods <- c("BY", "BH", "OrderShapeEM", "SABHA")

data.obj <- list()
res <- list()

# n.rep * length(alphas)
for (i in 1: length(methods)){
  res[[methods[i]]]$fdp <- matrix(NA, n.rep, length(alphas))
  res[[methods[i]]]$pd <- matrix(NA, n.rep, length(alphas))
}

for (i in 1: n.rep){
  cat("Replicate ", i, ":\n")
  data.obj[[i]] <- SimuData(tau       = tau,
                            ipt       = ipt,    
                            n.gene    = n.gene,
                            sig.str   = sig.str,
                            n.signal  = n.signal,
                            n.sample  = n.sample,
                            info.str  = info.str)
  pvalue <- data.obj[[i]]$pvalue
  x      <- data.obj[[i]]$x
  truth  <- data.obj[[i]]$truth
  if (length(x)<n.gene)
    x <- c(x, rep(1, (n.gene-length(x))))
  x[is.na(x)] <- 1
  
  # hist(pvalue, main = "Primary p-values", xlab = "p-values")
  
  # BY (original SPARK)
  pval.by <- p.adjust(pvalue, method = "BY")
  
  # BH
  pval.bh <- p.adjust(pvalue, method = "BH")
  
  # OrderShapeEM
  res.OrderShapeEM  <- OrderShapeEM(pvalue, x, OrderShapeEM.control(trace = TRUE))
  
  # SABHA
  tau.sabha = 0.5; 
  eps = 0.1 # parameters for SABHA
  index0 <- order(x)
  index1 <- order(index0)
  pvals  <- pvalue[index0]
  qhat   <- Solve_q_step(pvals, tau.sabha, eps) # compute qhat for SABHA
  
  for (j in 1:length(alphas)){
    alpha <- alphas[j]
    
    # BY procedure
    res$BY$fdp[i,j] <- sum(pval.by <= alpha & !truth)/max(sum(pval.by <= alpha), 1)
    res$BY$pd[i,j]  <- sum(pval.by <= alpha & truth) / sum(truth)
    
    # BH procedure
    res$BH$fdp[i,j] <- sum(pval.bh <= alpha & !truth)/max(sum(pval.bh <= alpha), 1)
    res$BH$pd[i,j]  <- sum(pval.bh <= alpha & truth) / sum(truth)
    
    # OrderShapeEM
    res$OrderShapeEM$fdp[i,j] <- sum(res.OrderShapeEM$fdr <= alpha & !truth)/max(sum(res.OrderShapeEM$fdr <= alpha), 1)
    res$OrderShapeEM$pd[i,j]  <- sum(res.OrderShapeEM$fdr <= alpha & truth) / sum(truth)
    
    # SABHA
    rej.sabha      <- SABHA_method(pvals, qhat, alpha, tau.sabha)
    rej.sabha.vec  <- rep(0, n.gene)
    rej.sabha.vec[rej.sabha] <- 1
    rej.sabha.vec <- rej.sabha.vec[index1]
    res$SABHA$fdp[i,j] <- sum(rej.sabha.vec & !truth) / max(sum(rej.sabha.vec), 1)
    res$SABHA$pd[i,j]  <- sum(rej.sabha.vec & truth) / sum(truth)
  }
}

for (k in 1:length(methods)){
  res[[k]]$fdr <- colMeans(res[[k]]$fdp)
  res[[k]]$pwr <- colMeans(res[[k]]$pd)
}


data.name <- paste0("./simulation/",info.str, "_Sig", n.signal/n.gene,"_N",n.gene,".RData")
save(data.obj,res, file = data.name)