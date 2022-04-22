
RunAll <- function(pval, covar, alpha){
  # run all methods: BY, BH, ST, SABHA, OrderShapeEM
  source('./R/All_q_est_functions.R')
  require(qvalue)
  require(Iso)
  source('./R/OrderShapeEM.R')
  
  # BY (original SPARK)
  pval.by <- p.adjust(pval, method = "BY")
  rej.by <- which(pval.by <= alpha)
  
  # BH
  pval.bh <- p.adjust(pval, method = "BH")
  rej.bh <- which(pval.bh <= alpha)
  
  # OrderShapeEM
  orderfdr.obj <- OrderShapeEM(pval, covar, OrderShapeEM.control(trace = FALSE))
  sum(orderfdr.obj$fdr < alpha)
  # 1112 SE identified by OrderShapeEM
  rej.osem <- which(orderfdr.obj$fdr <= alpha)
  
  # SABHA
  tau.sabha = 0.5; 
  eps = 0.1 # parameters for SABHA
  pvals  <- pval[order(covar)]
  qhat   <- Solve_q_step(pvals, tau.sabha, eps)
  rej.sabha <- SABHA_method(pval, qhat, alpha, tau.sabha)
  
  return(Rej.res = list(rej.by = rej.by, rej.bh = rej.bh, rej.sabha = rej.sabha, rej.osem = rej.osem))
}