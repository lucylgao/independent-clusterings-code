library(methods)
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

# makes Pi according to equation 4.14
make_Pi_vary <- function(delta, K) { 
  delta * diag(K) / K + (1 - delta) * matrix(1, K, K) / K^2 
}


##### Functions for running the simulation study ###### 
counter <- 0
# One dataset in the simulation study 
do_one <- function(delta, n) { 
  Pi <- make_Pi_vary(delta, 4) 
  
  mu1 <- cbind(c(2, -2), c(2, -1), c(-2, 1), c(-2, 2))
  mu2 <- cbind(c(-2, -2), c(-2, -1), c(2, 1), c(2, 2))
  
  # Simulate dataset
  dat <- generate_Gaussian_data(n, Pi, mu1, mu2, 0.4)
  x <- dat$x 
  cl <- dat$cl
  
  results <- rep(0, 6) 
  invisible(capture.output(results[1] <- tryCatch(test_indep_clust(x, K1=2, K2=2,  
                                                                   model1="EII", model2="EII", init1="EII", init2="EII")$pval, 
                                                  error=function(cond) { return(NA) } ) ))
  invisible(capture.output(results[2] <- tryCatch(test_indep_clust(x, K1=4, K2=4,  
                                                                   model1="EII", model2="EII", init1="EII", init2="EII")$pval, 
                                                  error=function(cond) { return(NA) } ) ))
  invisible(capture.output(results[3] <- tryCatch(test_indep_clust(x, K1=6, K2=6,  
                                                                   model1="EII", model2="EII", init1="EII", init2="EII")$pval, 
                                                  error=function(cond) { return(NA) } ) )) 
  invisible(capture.output(use.bic <- tryCatch(test_indep_clust(x, 
                                                                model1="EII", model2="EII", init1="EII", init2="EII"), 
                                               error=function(cond) { return(NA) } ) ))
  results[4] <- use.bic$pval
  results[5] <- use.bic$K1 
  results[6] <- use.bic$K2
  
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
rownames(results) <- c("K2-pval", "K4-pval", 
                       "K6-pval", "KBIC-pval", 
                       "KBIC-V1", "KBIC-V2")
save(results, file=paste("chooseK-setting1-n", n, "-delta", delta, 
                         "-seed", seed, ".Rdata", sep=""))