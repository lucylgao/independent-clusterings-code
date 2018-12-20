##### Functions for simulating multi-view data ###### 
# Generates n multi-view cluster memberships given a Pi matrix
generate_memberships <- function(n, pi) {
  stopifnot(pi >= 0, abs(sum(pi) - 1) < 1e-15)
  ii <- sample(length(pi), n, replace = TRUE, prob = pi)
  k1 <- nrow(pi)
  k2 <- ncol(pi)
  cpairs <- cbind(rep(1:k1, k2), rep(1:k2, each = k1))
  cpairs[ii, ]
}

# Given multi-view cluster memberships, 
# multi-view cluster means (mu1 & mu2) 
# and a shared sigma for the sig^2*I covariance structure, 
# generates multi-view Gaussian mixture data of sample size n
generate_data_from_memberships <- function(cl, mu1, mu2, Sig) {
  n <- nrow(cl)
  p1 <- nrow(mu1)
  p2 <- nrow(mu2)
  list(matrix(rnorm(n * p1), n, p1)%*%Sig + t(mu1)[cl[, 1], ],
       matrix(rnorm(n * p2), n, p2)%*%Sig + t(mu2)[cl[, 2], ])
}

# Given multi-view cluster means (mu1 & mu2), 
# a shared sigma for the sig^2*I covariance structure, 
# and a probability matrix Pi, 
# generates multi-view Gaussian mixture data & memberships 
# of sample size n 
generate_Gaussian_data = function(n, pi, mu1, mu2, Sig) { 
  cl <- generate_memberships(n, pi)
  x <- generate_data_from_memberships(cl, mu1, mu2, Sig)
  return(list(x=x, cl=cl))
}

make_Pi_vary <- function(delta, K) { 
  delta * diag(K) / K + (1 - delta) * matrix(1, K, K) / K^2 
}

meila_mclust <- function(x, K1 = NULL, K2 = NULL, model1, model2,
                         init1, init2, B=200) { 
  n <- dim(x[[1]])[1]
  
  if(is.null(K1)) {
    EM.View1 <- mclust::Mclust(x[[1]], G=2:9, modelNames=c(model1),
                               initialization=list(hcPairs=hc(x[[1]],
                                                              modelName=init1)))
    K1 <- EM.View1$G
  } else { 
    EM.View1 <- mclust::Mclust(x[[1]], K1, modelNames=c(model1), 
                             initialization=list(hcPairs=hc(x[[1]], 
                                                            modelName=c(init1))))
  }
  
  if(is.null(K2)) {
    EM.View2 <- mclust::Mclust(x[[2]], G=2:9, modelNames=c(model2),
                               initialization=list(hcPairs=hc(x[[2]],
                                                              modelName=init2)))
    K2 <- EM.View2$G
  } else { 
    EM.View2 <- mclust::Mclust(x[[2]], K2, modelNames=c(model2), 
                             initialization=list(hcPairs=hc(x[[2]], 
                                                            modelName=c(init2))))
  } 
  
  # Compute mutual information test statistic
  clustTable <- table(EM.View1$classification, EM.View2$classification)
  piEst <- clustTable/n 
  piEst1 <- rowSums(clustTable)/n 
  piEst2 <- colSums(clustTable)/n
  
  meilaStat <- 2*sum(clustTable*log(t(t(piEst/piEst1)/piEst2)), na.rm=TRUE)
  
  # Obtain p-value with permutation approach 
  meilaStatp <- rep(0, B)
  for(b in 1:B) { 
    clustTablep <- table(EM.View1$classification, EM.View2$classification[sample(1:n, replace=F)])
    piEstp <- clustTablep/n 
    
    meilaStatp[b] <- 2*sum(clustTablep*log(t(t(piEstp/piEst1)/piEst2)), na.rm=TRUE)
  }
  
  # Obtain p-value with chi^2 approximation
  pvalChi <- 1 - pchisq(meilaStat, df=(K1-1)*(K2-1))
  pvalp <- mean(ifelse(meilaStatp >= meilaStat, 1, 0))
  
  # Return test statistic and p-values
  return(c(meilaStat, pvalChi, pvalp))
}

##### Functions for running the simulation study ###### 
counter <- 0
# One dataset in the simulation study 
do_one <- function(delta, n) { 
  Pi <- make_Pi_vary(delta, 3) 
  
  mu1 <- cbind(c(0, 2), c(0, -2), c(sqrt(12), 0))
  mu2 <- cbind(c(-2, 0), c(0, sqrt(12)), c(2, 0))
  
  # Var mat = 2.25 & 0.5 // 0.5 & 2.25
  Sig <- chol(rbind(c(1.5^2, 0.5), c(0.5, 1.5^2)))
  
  # Simulate dataset
  dat <- generate_Gaussian_data(n, Pi, mu1, mu2, Sig)
  x <- dat$x 
  cl <- dat$cl
  
  results <- rep(NA, 4)
  invisible(capture.output(results[1] <- tryCatch(test_indep_clust(x, K1=3, K2=3,  
                                                                   model1="EEE", model2="EEE", init1="EEE", init2="EEE")$pval, 
                                                  error=function(cond) { return(NA) } ) ))
  invisible(capture.output(results[2] <- tryCatch(test_indep_clust(x,   
                model1="EEE", model2="EEE", init1="EEE", init2="EEE")$pval, 
                error=function(cond) { return(NA) } ) ))
  
  invisible(capture.output(results[3] <- tryCatch(meila_mclust(x, 
                                                  K1=3, K2=3, model1="EEE", model2="EEE", init1="EEE", init2="EEE")[3], 
                                                  error=function(cond) { return(NA) } ) ))
  
  invisible(capture.output(results[4] <- tryCatch(meila_mclust(x, model1="EEE", model2="EEE", init1="EEE", init2="EEE")[3], 
                                                  error=function(cond) { return(NA) } ) ))
  
  # Print progress updates once ~25%, 50%, 75%, and 100% of the job is done
  counter <- counter + 1; 
  if(counter %in% c(42, 84, 125, 167)) 
    cat(paste(round((counter/167)*100), "% has been processed"), sep="")
  
  # Return results
  return(results)
}


##### Running simulation study ###### 
# iterate or parallelize
# through the following settings: 
# delta \in seq(0, 1, length=9), n \in {25, 50, 100, 250, 500}, seed \in 1:12
set.seed(seed, kind = "L'Ecuyer-CMRG")
results <- replicate(167, do_one(delta, n))
save(results, file=paste("dense-gaussian-results-n", n, "-delta", delta, 
                         "-seed", seed, ".Rdata", sep="")) 

