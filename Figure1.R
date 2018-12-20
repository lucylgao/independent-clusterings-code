###### data generation functions ####### 
#' Generate n cluster membership pairs given a pi matrix
#' 
#' Generates an n-by-2 matrix of cluster membership pairs, (k, k'), where k is the
#' membership according to the first mode and k' is the membership according to the 
#' second mode.
#' 
#' @param n sample size
#' @param pi k1-by-k2 matrix of probabilities
generate_memberships <- function(n, pi) {
  stopifnot(pi >= 0, sum(pi) == 1)
  ii <- sample(length(pi), n, replace = TRUE, prob = pi)
  k1 <- nrow(pi)
  k2 <- ncol(pi)
  cpairs <- cbind(rep(1:k1, k2), rep(1:k2, each = k1))
  cpairs[ii, ]
}

#' Generate data based on cluster memberships
#' @param cl n-by-2 matrix of cluster memberships
#' @param mu1 p1-by-k1 matrix of mode-1 mean vectors
#' @param mu2 p2-by-k2 matrix of mode-2 mean vectors
#' @param sig standard deviation of Gaussian (a scalar and shared 
#'        across all clusters)
generate_data_from_memberships <- function(cl, mu1, mu2, sig = 0.1) {
  n <- nrow(cl)
  p1 <- nrow(mu1)
  p2 <- nrow(mu2)
  list(matrix(sig * rnorm(n * p1), n, p1) + t(mu1)[cl[, 1], ],
       matrix(sig * rnorm(n * p2), n, p2) + t(mu2)[cl[, 2], ])
}


###### make figure ####### 
#' Make row of plots for a given pi
#' 
#' When p1 = p2 = 2.
#' 
#' Colors show memberships according to first view; shapes show memberships
#' according to second view
#' 
#' @param pi a k1 by k2 matrix of probabilities
#' @param mu1 k1-by-2 matrix
#' @param mu2 k2-by-2 matrix
new_make_row <- function(pi, mu1, mu2) {
  # generate some data to show:
  cl <- generate_memberships(n = 15, pi = pi)
  x <- generate_data_from_memberships(cl, mu1, mu2, sig = 0.3)
  k1 <- nrow(pi); k2 <- ncol(pi)
  # do plotting:
  par(mfrow = c(1, 2), mar = rep(0.5, 4))
  if (k1 == 2)
    cols <- c("#525252", "#cccccc")
  if (k1 == 3)
    cols <- c("blue", "orange", "green")
  for (i in 1:2) {
    plot(x[[i]], xlab = "", ylab = "", 
         pch = c(17, 19)[cl[, 2]], cex = 4, col = cols[cl[, 1]], axes = FALSE,
         xlim = c(-0.7, 2.7), ylim = c(-0.7, 2.7))
    points(x[[i]], pch = c(17, 19)[cl[, 2]], cex = 4.5, col = 1) # outline
    points(x[[i]], pch = c(17, 19)[cl[, 2]], cex = 4, col = cols[cl[, 1]])
    box()
  }
}

# rank 1 pi: Figure 1(i)
pi <- matrix(1, 2, 2) / 4
mu1 <- cbind(c(2, 2), c(0.3, 0))
mu2 <- cbind(c(1, 0), c(1, 2))
set.seed(122)
png("Figure1(i).png", width = 8, height = 4, units = 'in', res = 300)
new_make_row(pi, mu1, mu2)
dev.off()

# diagonal pi: Figure 1(ii)
pi <- diag(2) / 2
mu1 <- cbind(c(2, 2), c(0.3, 0))
mu2 <- cbind(c(1, 0), c(1, 2))
set.seed(122)
png("Figure1(ii).png", width = 8, height = 4, units = 'in', res = 300)
new_make_row(pi, mu1, mu2)
dev.off()


# intermediate pi: Figure 1(iii)
pi <- (matrix(1, 2, 2) / 4)*0.5 + diag(2)*0.25
mu1 <- cbind(c(2, 2), c(0.3, 0))
mu2 <- cbind(c(1, 0), c(1, 2))
set.seed(122)
png("Figure1(iii).png", width = 8, height = 4, units = 'in', res = 300)
new_make_row(pi, mu1, mu2)
dev.off()


effective_rank <- function(M) { 
  M.sv <- svd(M)$d
  return(sum(M.sv)/max(M.sv))
}

effective_rank(pi)


