
pdf("FigureS2(i).pdf", width=8, height=4)
par(mfrow=c(1, 2))
mu1 <- cbind(c(2, -2), c(2, -1), c(-2, 1), c(-2, 2))
mu2 <- cbind(c(-2, -2), c(-2, -1), c(2, 1), c(2, 2))
plot(mu1[1, ], mu1[2, ], 
     pch=as.character(1:4), xlab="", ylab="", 
     main="View 1", cex=1.5, 
     xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5))
plot(mu2[1, ], mu2[2, ], 
     pch=as.character(1:6), xlab="", ylab="", 
     main="View 2", 
     cex=1.5, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5))
dev.off()

pdf("FigureS2(ii).pdf", width=8, height=4)
par(mfrow=c(1, 2))
mu2 <- mu2[, c(4, 1, 2, 3)]
plot(mu1[1, ], mu1[2, ], 
     pch=as.character(1:4), xlab="", ylab="", 
     main="View 1", cex=1.5, 
     xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5))
plot(mu2[1, ], mu2[2, ], 
     pch=as.character(1:6), xlab="", ylab="", 
     main="View 2", 
     cex=1.5, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5))
dev.off()
