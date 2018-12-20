# Run sim-study-sec4-sec5.R first to generate simulation results for specific
# simulation set-ups.
# The following file collates the simulation results corresponding 
# to those set-ups in an easy to read format.

# For Figures 2&3, set: p = 10, K = 6
# For Figures S3&S6, set: p = 10, K = 3
# For Figures S4&S7, set: p = 100, K = 3
# For Figures S5&S8, set: p = 100, K = 6

# Constructs and saves a data frame containing all of the knots 
# (delta ranging from 0 to 1, n ranging from 25 to 500) in the 
# three power curves corresponding to different sig values
# in the simulation setup 
sig.range <- if(p == 10) c(2.4, 4.8, 9.6) else c(4.8, 9.6, 19.2)
offset <- ifelse(K == 6, 3,1)


delta_seq <- seq(0, 1, length.out=9)
power_cat <- c()

nSim <- 12*167
for(sig in sig.range) { 
  for(n in c(25, 50, 100, 250, 500)) { 
    for(delta_ind in 1:9) { 
      delta <- delta_seq[delta_ind]
      results_cat <- c()
      for(seed in 1:12) { 
        load(file=paste("results-K", K, "-p", p, "-sig", sig, "-n", n, 
                        "-delta", delta, "-seed", seed, ".Rdata", sep=""))
        results_cat <- cbind(results_cat, results)
      }
      
            power_cat <- rbind(power_cat, rbind(
        c("test_indep_per", mean(ifelse(results_cat[2, ] < 0.05, 1, 0)), 
          sqrt(mean(ifelse(results_cat[2, ] < 0.05, 1, 0))*
                 (1 - mean(ifelse(results_cat[2, ] < 0.05, 1, 0)))/nSim), 
          sig, n, K, K, delta, p), 
        c(paste("test_indep_perK", K-offset, sep=""), 
          mean(ifelse(results_cat[4, ] < 0.05, 1, 0),  na.rm=T), 
          sqrt(mean(ifelse(results_cat[4, ] < 0.05, 1, 0), 
                    na.rm=T)*
                 (1 - mean(ifelse(results_cat[4, ] < 0.05, 1, 0), 
                           na.rm=T))/(nSim- sum(is.na(results_cat[4, ])))), 
          sig, n, K, K, delta, p),
        c(paste("test_indep_perK", K+offset, sep=""), 
          mean(ifelse(results_cat[6, ] < 0.05, 1, 0), na.rm=T), 
          sqrt(mean(ifelse(results_cat[6, ] < 0.05, 1, 0),  na.rm=T)*
                 (1 - mean(ifelse(results_cat[6, ] < 0.05, 1, 0), 
                           na.rm=T))/(nSim- sum(is.na(results_cat[6, ])))), 
          sig, n, K, K, delta, p),
        c("meila_mclust", mean(ifelse(results_cat[8, ] < 0.05, 1, 0)), 
          sqrt(mean(ifelse(results_cat[8, ] < 0.05, 1, 0))*
                 (1 - mean(ifelse(results_cat[8, ] < 0.05, 1, 0)))/nSim), 
          sig, n, K, K, delta, p),
        c("meila_mclust_per", mean(ifelse(results_cat[9, ] < 0.05, 1, 0)), 
          sqrt(mean(ifelse(results_cat[9, ] < 0.05, 1, 0))*
                 (1 - mean(ifelse(results_cat[9, ] < 0.05, 1, 0)))/nSim), 
          sig, n, K, K, delta, p),
        c("ARI", mean(ifelse(results_cat[11, ] < 0.05, 1, 0)), 
          sqrt(mean(ifelse(results_cat[11, ] < 0.05, 1, 0))*
                 (1 - mean(ifelse(results_cat[11, ] < 0.05, 1, 0)))/nSim), 
          sig, n, K, K, delta, p)))
    }
  }
  power_cat <- data.frame(power_cat)
  names(power_cat) <- c("Method", "power", "se", "sig", "n", 
                        "K1", "K2", "delta", "p")
  power_cat$power <- as.numeric(as.character(power_cat$power))
  power_cat$n <- as.numeric(as.character(power_cat$n))
  power_cat$se <- as.numeric(as.character(power_cat$se))
  power_cat$sig <- as.numeric(as.character(power_cat$sig))
  power_cat$K1 <- as.numeric(as.character(power_cat$K1))
  power_cat$K2 <- as.numeric(as.character(power_cat$K2))
  power_cat$delta <- as.numeric(as.character(power_cat$delta))
  power_cat$p <- as.numeric(as.character(power_cat$p))
  
  df <- power_cat
  cat(paste("sig=",sig, " NA= ", sum(is.na(power_cat))), sep="")
  save(df, file=paste("df-K", K, "-p", p, "-sig", sig, ".Rdata", sep=""))
  power_cat <- c()
}



