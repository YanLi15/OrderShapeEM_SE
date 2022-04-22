######################################################################
#               Real data guided Simulation function                 #
#  SPARK primary p-values + two-sample RNA-seq auxiliary covariates  #
######################################################################
# @param tau: a value in (0,1), the variance of non-spatial residual error. median 0.35 as default.
# @param ipt: a value in {2, 3, 4}, the index of pattern. 2 as default.
# @param n.gene: a positive integer, the total number of simulated genes.
# @param sig.str: a string, strength of SE gene signals. Moderate as default, also support Weak and Strong.
# @param n.signal: a positive interger no larger than n.gene, the number of SE genes.
# @param n.sample: a positive interger, the number of samples under each condition in RNA-seq data.
# @param info.str: a string, auxiliary informativeness. Uninformative, Weak, or Strong.

SimuData <-
  function(tau      = 0.35,
           ipt      = 2,
           n.gene   = 10000,
           sig.str  = "Moderate",
           n.signal = 1000,
           n.sample = 50,
           info.str = c("Uninformative", "Weak", "Strong")) {
    require(SPARK)
    require(compcodeR)
    require(DESeq2)
    
    ## Generate the spatial count data based on the SPARK results from the real
    ## mouse olfactory bulb data
    load("./Output/MOB.RData")
    
    beta <- sapply(1:length(spark_MOB@res_vc), function(x) {
      spark_MOB@res_vc[[x]]$coefficients
    })
    nb <- sapply(1:4, function(x) {
      log(x * exp(median(beta)))
    })
    load("./RealData/MOB_Pattern_SpatialDE.rds")
    info  <-
      read.csv("./RealData/Rep11_MOB_info_spark.csv",
               row.names = 1)
    newN  <- info$total_counts
    pattern <- datlist[[paste0("pattern", ipt)]]
    grp <- as.numeric(pattern[, 2] > mean(pattern[, 2])) + 1
    
    if (sig.str == "Weak")
      ifc <- 2
    if (sig.str == "Moderate")
      ifc <- 3
    if (sig.str == "Strong")
      ifc <- 4
    uu  <- c(nb[1], nb[ifc], nb[1])
    
    truth <- c(rep(1, n.signal), rep(0, (n.gene - n.signal)))
    lambda <- sapply(1:n.gene, function(x) {
      truth[x] * exp(uu[grp] + rnorm(length(uu[grp]), 0, tau)) + (!truth[x]) *
        exp(uu[3] + rnorm(length(uu[grp]), 0, tau))
    })
    newCt <- lapply(1:(n.gene), function(x) {
      rpois(length(lambda[, x]), newN * lambda[, x])
    })
    countdata <- data.frame(do.call(rbind, newCt))
    rownames(countdata) <- paste0("gene", 1:nrow(countdata))
    colnames(countdata) <- pattern[, 1]
    
    ## Compute primary p-values with SPARK
    spark_MOB <-
      CreateSPARKObject(
        counts = countdata,
        location = info[, 1:2],
        percentage = 0.1,
        min_total_counts = 10
      )
    
    spark_MOB@lib_size <- info$total_counts
    spark_MOB <-
      spark.vc(
        spark_MOB,
        covariates = NULL,
        lib_size = spark_MOB@lib_size,
        num_core = 1,
        verbose = T,
        fit.maxiter = 500
      )
    spark_MOB <- spark.test(spark_MOB, check_positive = T, verbose = T)
    pvalue <- spark_MOB@res_mtest[, "combined_pvalue"]
    
    closeAllConnections()
    
    ## Generate two-sample RNA-seq data and compute auxiliary covariates
    RNAseq <-
      generateSyntheticData(
        dataset = "RNASeq",
        n.vars  = n.gene,
        samples.per.cond = n.sample,
        n.diffexp = n.signal,
        seqdepth = 1e7,
        fraction.upregulated = 0.5,
        filter.threshold.total = 1,
        effect.size = 1.5)
    
    countdata <- RNAseq@count.matrix
    group <-
      factor(c(rep("condition1", n.sample), rep("condition2", n.sample)), levels = c("condition1", "condition2"))
    colData <-
      data.frame(row.names = colnames(countdata), group = group)
    dds <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData   = colData,
                                  design    = ~ group)
    # dispersion estimates and differential analysis
    dds <- DESeq(dds)
    res <- results(dds)
    aux.pvalue <- res@listData$pvalue
    aux.truth <- RNAseq@variable.annotations$differential.expression
    
    if (info.str == "Uninformative") {
      aux.pvalue <- runif(n.gene, 0, 1)
      aux.truth <- rep(0, n.gene)
    }
    if (info.str == "Weak") {
      index <- sample(1:n.gene)
      
      aux.pvalue <- aux.pvalue[index]
      aux.truth  <- aux.truth[index]
    }
    if (info.str == "Medium"){
      if (n.signal * 5 < n.gene)
        index <-
          c(sample(1:(n.signal * 5)), (n.signal * 5 + 1):n.gene)
      else
        index <- sample(1:n.gene)
      
      aux.pvalue <- aux.pvalue[index]
      aux.truth  <- aux.truth[index]
    }
    if (info.str == "Strong") {
      if (n.signal * 2 < n.gene)
        index <-
          c(sample(1:(n.signal * 2)), (n.signal * 2 + 1):n.gene)
      else
        index <- sample(1:n.gene)
      
      aux.pvalue <- aux.pvalue[index]
      aux.truth  <- aux.truth[index]
    }
    
    return(list(
      pvalue = pvalue,
      x = aux.pvalue,
      truth = truth,
      aux.truth = aux.truth
    ))
  }
