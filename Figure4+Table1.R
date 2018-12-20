library(impute)
library(multiviewtest)

####### Reading in the raw data #######
# participant demographics 
participants <- read.csv("participant_data.csv")[, -c(1)]

# clinical data
clinical <- read.table("clinical_labs_DS.txt", sep="\t", 
                       comment.char="", quote="",
                       stringsAsFactors=F,  header=T)
# pick out the binary variables
clinical <- clinical[!clinical$VARNAME %in% 
                       unique(clinical[clinical$UNIT
                                       %in% c("bool", "pattern"), "VARNAME"]), ]
clinical <- clinical[, c(1, 2, 3, 9)]
names(clinical) <- c("username", "time", "varname", "value")
clinical <- reshape(clinical,idvar=c("username", "time"),  
                    timevar=c("varname"), direction="wide")
clinical <- merge(clinical, participants)[, c(1:2, 210:212, 3:209)]
clinical$username <- as.numeric(matrix(unlist(strsplit(clinical$username, "P")), 
                                       2, nrow(clinical))[2, ])
clinical <- clinical[order(clinical$username, clinical$time), ]
names(clinical)[1] <- "id"
rownames(clinical) <- NULL

# protein data
protein <- read.table("proteomics_DS.txt", sep="\t", 
                      comment.char="", quote="",
                      stringsAsFactors=F, header=T)
protein$PEPTIDE[protein$PEPTIDE == ""] <- "None"
protein <- protein[order(protein$VARNAME,protein$CATEGORY, 
                         protein$PEPTIDE), ]
protein$VARNAME <- paste(protein$VARNAME, 
                         protein$CATEGORY, 
                         protein$PEPTIDE, sep=".")
protein <- protein[, c(1, 2, 3, 11)]
names(protein) <- c("username", "time", "varname", "value")
protein <- reshape(protein,idvar=c("username", "time"),  
                   timevar=c("varname"),  direction="wide")
protein <- merge(protein, participants)[, c(1:2, 271:273, 3:270)]
protein$username <- as.numeric(matrix(unlist(strsplit(protein$username, "P")), 
                                      2, nrow(protein))[2, ])
protein <- protein[order(protein$username, protein$time), ]
names(protein)[1] <- "id"
rownames(protein) <- NULL

# metabolite data
metab <- read.table("metabolites_DS.txt", sep="\t", 
                    comment.char="", quote="",
                    stringsAsFactors=F,  header=T)
metab <- metab[, c(1, 2, 3, 11)]
names(metab) <- c("username", "time", "varname", "value")
metab <- reshape(metab,idvar=c("username", "time"),  
                 timevar=c("varname"),   direction="wide")
metab <- merge(metab, participants)[, c(1:2, 645:647, 3:644)]
metab$username <- as.numeric(matrix(unlist(strsplit(metab$username, "P")), 
                                    2, nrow(metab))[2, ])
metab <- metab[order(metab$username, metab$time), ]
names(metab)[1] <- "id"
rownames(metab) <- NULL

####### Function for pre-processing the data ####### 
preprocess_data <- function(view, miss) { 
  cat(paste("Originally: ", nrow(view), " subjects, ", 
            ncol(view) - 5, " features.\n", sep=""))
  
  # Only keep features which exist in 75% or more subjects
  miss.subj <- apply(view[, -c(1:5)], 2, function(x) sum(is.na(x)))/nrow(view)
  cat(paste("Removed ", sum(!miss.subj < miss), 
            " features for being available in < ", 
            (1 - miss)*100, "% subjects.\n", sep=""))  
  view <- view[, c(1:5, which(miss.subj < miss)+5)]
  
  # Only keep subjects that have 75% or more features
  miss.feat <- apply(view[, -c(1:5)], 1, function(x) sum(is.na(x)))/
    ncol(view[, -c(1:5)])
  cat(paste("Removed ", sum(!miss.feat < miss), " out of ", 
            nrow(view), " subjects for having < ", (1 - miss)*100, 
            "% features available.\n", sep=""))
  view <- view[miss.feat < miss, ]
  
  # Only keep features with more than 0 sd 
  sd0 <- apply(view[, -c(1:5)], 2, function(x) sd(x, na.rm=T)) == 0
  cat(paste("Removed", sum(sd0), "features for having 0 sd.\n" ))
  view <- view[, c(1:5, which(!sd0)+5)]
  
  # Imputing the missing values 
  cat(paste("Imputing", sum(is.na(view[, -c(1:5)])), "values out of", 
            nrow(view)*ncol(view[, -c(1:5)]), "values ... \n"))
  view[, -c(1:5)] <- impute.knn(as.matrix(view[, -c(1:5)]))$data
  
  # Regressing out gender
  cat("Regressing out gender ... \n")
  for(i in 6:ncol(view)) { 
    if(!all(view[!is.na(view[, i]), "gender"] == 
            view[!is.na(view[, i]), "gender"][1])) 
      view[, i] <- resid(lm(view[, i]~ifelse(view$gender == "F", 1, 0), 
                            na.action=na.exclude))
  }
  
  # Scale the view
  cat(paste("Scaling the view ... \n"))
  view[, -c(1:5)] <- scale(view[, -c(1:5)])
  
  cat(paste("Final: ", nrow(view), " subjects, ", 
            ncol(view) - 5, " features.\n", sep=""))
  
  return(view)
}

###### P100 Analysis ###### 
# Function to compare clusterings and plot/print the results
comparing_clusters <- function(view1, view2, K1=NULL, K2=NULL, 
                               model="EII", miss, B=10^4, step=0.001,
                               main1, main2, titleplot) { 
  view1 <- preprocess_data(view1, miss)
  view2 <- preprocess_data(view2, miss)
  
  
  # Cluster full view 1 and full view 2 
  # Compare clusterings on the subsets of view 1 and view 2 such that 
  # observations exist in both view 1 and view 2
  subset1 <- which(view1$id %in% intersect(view1$id, view2$id))
  subset2 <- which(view2$id %in% intersect(view1$id, view2$id))
  results <- test_indep_clust_subset(x=list(as.matrix(view1[, -c(1:5)]),
                               as.matrix(view2[, -c(1:5)])), K1=K1, K2=K2, 
                          model1=model, model2=model, B=B, 
                          subset1=subset1, 
                          subset2=subset2, 
                          init1="EII", init2="EII")
  
  # PCA Plots
  pca1 <- prcomp(view1[, -c(1:5)])
  pca2 <- prcomp(view2[, -c(1:5)])
  
  cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
  pchs <- c(17, 19, 15, 18, 25, 8, 3)
  pdf(paste(titleplot, ".pdf", sep=""),  width=12, height=12)
  par(mfrow=c(1, 2), pty="s")
  v1.ids <- view1$id %in% intersect(view1$id, view2$id)
  v2.ids <- view2$id %in% intersect(view1$id, view2$id)
  plot(pca1$x[subset1, ], cex=2,  col=cols[results$modelfit1$classification[subset1]], 
       pch=pchs[results$modelfit2$classification[subset2]], 
       bg=cols[results$modelfit1$classification[subset1]], 
       axes=FALSE, xlab="", ylab="")
  title(xlab=paste("Principal Component 1 (", 
                   round(summary(pca1)$importance[2, 1]*100, 2), "%)", sep=""), 
        ylab=paste("Principal Component 2 (", 
                   round(summary(pca1)$importance[2, 2]*100, 2), "%)", sep=""), line=0.5, 
        main=main1, cex.main=2, cex.lab=2)
  points(pca1$x[subset1, ], pch = pchs[results$modelfit2$classification[subset2]],
         cex = 2.4, col = "black",   bg="black") # outline
  points(pca1$x[subset1, ], pch = pchs[results$modelfit2$classification[subset2]], 
         cex = 2.3, col=cols[results$modelfit1$classification[subset1]], 
         bg=cols[results$modelfit1$classification[subset1]])
  box()
  plot(pca2$x[subset2, ],   cex=2, col=cols[results$modelfit1$classification[subset1]], 
       pch=pchs[results$modelfit2$classification[subset2]], 
       bg=cols[results$modelfit1$classification[subset1]],  
       axes=FALSE, xlab="", ylab="")
  title(xlab=paste("Principal Component 1 (", 
                   round(summary(pca2)$importance[2, 1]*100, 2), "%)", sep=""), 
        ylab=paste("Principal Component 2 (", 
                   round(summary(pca2)$importance[2, 2]*100, 2), "%)", sep=""), line=0.5, 
        main=main2, cex.main=2, cex.lab=2)
  points(pca2$x[subset2, ], pch=pchs[results$modelfit2$classification[subset2]],
         cex = 2.4, col = "black",  bg="black") # outline
  points(pca2$x[subset2, ], pch = pchs[results$modelfit2$classification[subset2]], 
         cex = 2.3, col=cols[results$modelfit1$classification[subset1]],
         bg=cols[results$modelfit1$classification[subset1]])
  box()  
  dev.off()
  
  return(list(PLRstat=results$PLRstat, Pi.est=results$Pi.est,
              K1=results$K1,  K2=results$K2, pval=results$pval,
              modelfit1=results$modelfit1, modelfit2=results$modelfit1, 
              n=length(subset1), p1=ncol(view1[-c(1:5)]), 
              p2=ncol(view2[-c(1:5)]),
              view1=view1, view2=view2))
}

# Generates Figure 4, and gets results for Table 1
set.seed(123)
# Generates Figure 4i
clin1clin3 <-comparing_clusters(clinical[clinical$time == 1, ], 
                                clinical[clinical$time == 3, ], 
                                miss=0.25,  main1="Clinical at Timepoint 1", 
                                main2="Clinical at Timepoint 3", 
                                titleplot="clin1clin3", B=10^5)
# Generates Figure 4ii
prot1prot3 <- comparing_clusters(protein[protein$time == 1, ], 
                                 protein[protein$time == 3, ], 
                                 miss=0.25,   main1="Proteomic at Timepoint 1", 
                                 main2="Proteomic at Timepoint 3", 
                                 titleplot="prot1prot3", B=10^5)
# Generates Figure 4iii
metab1metab3 <- comparing_clusters(metab[metab$time == 1, ], 
                                   metab[metab$time == 3, ], 
                                   miss=0.25, main1="Metabolomic at Timepoint 1", 
                                   main2="Metabolomic at Timepoint 3", 
                                   titleplot="metab1metab3", 
                                   K1=3, K2=3, B=10^5)

# Results for Table 1
clin1prot1 <- comparing_clusters(clinical[clinical$time == 1, ], 
                                 protein[protein$time == 1, ], 
                                 miss=0.25, main1="Clinical at Timepoint 1", 
                                 main2="Proteomic at Timepoint 1", 
                                 titleplot="clin1prot1", B=10^5)
clin1metab1 <- comparing_clusters(clinical[clinical$time == 1, ], 
                                  metab[metab$time == 1, ], 
                                  miss=0.25, K2 = 3, 
                                  main1="Clinical at Timepoint 1", 
                                  main2="Metabolomic at Timepoint 1", 
                                  titleplot="clin1metab1", B=10^5)
prot1metab1 <- comparing_clusters(protein[protein$time == 1, ], 
                                  metab[metab$time == 1, ], 
                                  miss=0.25, K2 = 3, 
                                  main1="Proteomic at Timepoint 1", 
                                  main2="Metabolomic at Timepoint 1", 
                                  titleplot="prot1metab1", B=10^5)

clin2prot2 <- comparing_clusters(clinical[clinical$time == 2, ], 
                                 protein[protein$time == 2, ], 
                                 miss=0.25, main1="Clinical at Timepoint 2", 
                                 main2="Proteomic at Timepoint 2", 
                                 titleplot="clin2prot2", B=10^5)
clin2metab2 <- comparing_clusters(clinical[clinical$time == 2, ], 
                                  metab[metab$time == 2, ], 
                                  miss=0.25, K2=4, main1="Clinical at Timepoint 2", 
                                  main2="Metabolomic at Timepoint 2", 
                                  titleplot="clin2metab2", B=10^5)
prot2metab2 <- comparing_clusters(protein[protein$time == 2, ], 
                                  metab[metab$time == 2, ], 
                                  miss=0.25, K2=3, main1="Proteomic at Timepoint 2", 
                                  main2="Metabolomic at Timepoint 2", 
                                  titleplot="prot2metab2", B=10^5)

clin3prot3 <- comparing_clusters(clinical[clinical$time == 3, ], 
                                 protein[protein$time == 3, ], 
                                 miss=0.25, main1="Clinical at Timepoint 3", 
                                 main2="Proteomic at Timepoint 3", 
                                 titleplot="clin3prot3", B=10^5)
clin3metab3 <- comparing_clusters(clinical[clinical$time == 3, ], 
                                  metab[metab$time == 3, ], 
                                  miss=0.25, main1="Clinical at Timepoint 3", 
                                  main2="Metabolomic at Timepoint 3", 
                                  titleplot="clin3metab3", 
                                  K2=3, B=10^5)
prot3metab3 <- comparing_clusters(protein[protein$time == 3, ], 
                                  metab[metab$time == 3, ], 
                                  miss=0.25, K2=3, main1="Proteomic at Timepoint 3", 
                                  main2="Metabolomic at Timepoint 3", 
                                  titleplot="prot3metab3", B=10^5)

