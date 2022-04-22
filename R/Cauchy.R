Cauchy <- function(Pvals){
    CCT <- colMeans(tan((0.5-Pvals)*pi))
    combined.pval <- 0.5-(atan(CCT))/pi
    return(combined.pval)
}