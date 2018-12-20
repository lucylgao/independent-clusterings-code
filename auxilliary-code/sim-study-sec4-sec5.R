library(multiviewtest)
library(mclust)

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
generate_data_from_memberships <- function(cl, mu1, mu2, sig) {
  n <- nrow(cl)
  p1 <- nrow(mu1)
  p2 <- nrow(mu2)
  list(matrix(sig * rnorm(n * p1), n, p1) + t(mu1)[cl[, 1], ],
       matrix(sig * rnorm(n * p2), n, p2) + t(mu2)[cl[, 2], ])
}

# Given multi-view cluster means (mu1 & mu2), 
# a shared sigma for the sig^2*I covariance structure, 
# and a probability matrix Pi, 
# generates multi-view Gaussian mixture data & memberships 
# of sample size n 
generate_Gaussian_data = function(n, pi, mu1, mu2, sig) { 
  cl <- generate_memberships(n, pi)
  x <- generate_data_from_memberships(cl, mu1, mu2, sig)
  return(list(x=x, cl=cl))
}

# Generate multi-view data for simulations 
# To generate data for each of the figures, iterate or parallelize
# through the following settings: 
# For all: delta \in seq(0, 1, length=9), n \in {25, 50, 100, 250, 500}, seed \in 1:12
# For Figure 2&3: p = 100, K = 6, sig \in {4.8, 9.6, 19.2}
# For Figure S1&S4: p = 10, K = 3, sig \in {2.4, 4.8, 9.6}
# For Figure S2&S5: p = 10, K = 6, sig \in {2.4, 4.8, 9.6}
# For Figure S3&S6: p = 100, K = 3, sig \in {4.8, 9.6, 19.2}
generate_multi_view_data <- function(delta, K, p, n, sig) { 
  Pi <- delta * diag(K) / K + (1 - delta) * matrix(1, K, K) / K^2 
  
  if(p == 2) { 
    mu1 <- cbind(c(2, 0), c(0, 2),  c(2, -2), c(-2, 0), c(0, -2), c(-2, 2))
    mu2 <- cbind(c(-2, 0), c(0, -2), c(-2, 2), c(2, 0), c(0, 2), c(2, -2)) 
  }
  
  if(p == 10) { 
    mu1 <- cbind(c(rep(2, 5), rep(0, 5)), 
                 c(rep(0, 5), rep(2, 5)), 
                 c(rep(2, 5), rep(-2, 5)), 
                 c(rep(-2, 5), rep(0, 5)), 
                 c(rep(0, 5), rep(-2, 5)), 
                 c(rep(-2, 5), rep(2, 5)))
    mu2 <- cbind(c(rep(-2, 6), rep(0, 4)), 
                 c(rep(0, 6), rep(-2, 4)), 
                 c(rep(-2, 4), rep(2, 6)),
                 c(rep(2, 6), rep(0, 4)), 
                 c(rep(0, 4), rep(2, 6)), 
                 c(rep(2, 4), rep(-2, 6)))
  } 
  
  if(p == 100) { 
    mu1 <- cbind(c(rep(2, 50), rep(0, 50)), 
                 c(rep(0, 50), rep(2, 50)), 
                 c(rep(2, 50), rep(-2, 50)), 
                 c(rep(-2, 50), rep(0, 50)), 
                 c(rep(0, 50), rep(-2, 50)), 
                 c(rep(-2, 50), rep(2, 50)))
    mu2 <- cbind(c(rep(-2, 60), rep(0, 40)), 
                 c(rep(0, 60), rep(-2, 40)), 
                 c(rep(-2, 40), rep(2, 60)),
                 c(rep(2, 60), rep(0, 40)), 
                 c(rep(0, 40), rep(2, 60)), 
                 c(rep(2, 40), rep(-2, 60)))
  } 
  
  
  mu1 <- mu1[, 1:K]
  mu2 <- mu2[, 1:K]  
  
  return(generate_Gaussian_data(n, Pi, mu1, mu2, sig))
}


##### Functions implementing alternative methods in Section 5 ###### 
# Implements the mutual information test described in Section 5 
meila_mclust <- function(x, K1, K2, model1, model2,
                         init1, init2, B=200) { 
  n <- dim(x[[1]])[1]
  
  # Model-based clustering of each view 
  EM.View1 <- mclust::Mclust(x[[1]], K1, modelNames=c(model1), 
                     initialization=list(hcPairs=hc(x[[1]], 
                                    modelName=c(init1))))
  EM.View2 <- mclust::Mclust(x[[2]], K2, modelNames=c(model2), 
                     initialization=list(hcPairs=hc(x[[2]], 
                                   modelName=c(init2))))
  
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

# Implements the ARI test described in Section 5 
ari <- function(x, K1, K2, model1, model2, init1, init2, B=200) { 
  n <- dim(x[[1]])[1]
  
  # Model-based clustering of each view 
  EM.View1 <- Mclust(x[[1]], K1, modelNames=c(model1), 
                     initialization=list(hcPairs=hc(x[[1]], 
                                        modelName=c(init1))))
  EM.View2 <- Mclust(x[[2]], K2, modelNames=c(model2), 
                     initialization=list(hcPairs=hc(x[[2]], 
                                        modelName=c(init2))))
  
  # Compute ARI 
  ari <- adjustedRandIndex(EM.View1$classification, EM.View2$classification)
  
  # Obtain p-value with permutation approach 
  arip <- replicate(B, adjustedRandIndex(EM.View1$classification, 
                      EM.View2$classification[sample(1:n, replace=F)]))
  
  pvalp <- mean(ifelse(arip >= ari, 1, 0))
  
  # Return test statistic and p-value
  return(c(ari, pvalp))
}

##### Functions for running the simulation study ###### 
counter <- 0
# One dataset in the simulation study 
do_one <- function(delta, K, p, n, sig) { 
  # If K = 3, also look at performance of K=2 and K=4 with PLRT 
  # Else K = 6, look at performance of K=3 and K=9 with PLRT
  offset <- ifelse(K == 3, 1, 3)
  
  # Simulate dataset
  dat <- generate_multi_view_data(delta, K, p, n, sig)
  x <- dat$x 
  cl <- dat$cl
  
  # Run PLRT under true K and misspecified K, mutual information test, 
  # and ARI test
  results <- rep(0, 11)
  results[1:2] <- tryCatch( 
    unlist(test_indep_clust(x, K1=K, K2=K,  
                            model1="EII", model2="EII", 
                            init1="EII", init2="EII")[c(1, 5)]), 
                            error=function(cond) { return(NA) } )
  results[3:4] <- tryCatch( 
    unlist(test_indep_clust(x, K1=K-offset, K2=K-offset,  
                            model1="EII", model2="EII", 
                            init1="EII", init2="EII")[c(1, 5)]), 
    error=function(cond) { return(NA) } )
  results[5:6] <-  tryCatch( 
    unlist(test_indep_clust(x, K1=K+offset, K2=K+offset,  
                            model1="EII", model2="EII", 
                            init1="EII", init2="EII")[c(1, 5)]), 
    error=function(cond) { return(NA) } )
  results[7:9] <- tryCatch(meila_mclust(x, K1=K, K2=K, 
                             model1="EII", model2="EII", 
                             init1="EII", init2="EII") , 
                           error=function(cond) { return(NA) } )    
  results[10:11] <- tryCatch(ari(x, K1=K, K2=K, 
                                 model1="EII", model2="EII", 
                                 init1="EII", init2="EII"), 
                             error=function(cond) { return(NA) } )    
  
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
# For all: delta \in seq(0, 1, length=9), n \in {25, 50, 100, 250, 500}, seed \in 1:12
# For Figures 2&3: p = 10, K = 6, sig \in {2.4, 4.8, 9.6}
# For Figures S3&S6: p = 10, K = 3, sig \in {2.4, 4.8, 9.6}
# For Figures S4&S7: p = 100, K = 3, sig \in {4.8, 9.6, 19.2}
# For Figures S5&S8: p = 100, K = 6, sig \in {4.8, 9.6, 19.2}
set.seed(seed, kind = "L'Ecuyer-CMRG")
results <- replicate(167, do_one(delta, K, p, n, sig))

# Save results
save(results, file=paste("results-K", K, "-p", p, "-sig", sig, "-n", n, 
                         "-delta", delta, "-seed", seed, ".Rdata", sep=""))
