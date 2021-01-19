# =================================
# = JGIL caller for allele counts =
# =================================

args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]

M <- matrix(c(4, 0, 0, 0, 3, 2, 1, 3, 2, 1, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 4, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 3, 2, 1, 3, 2, 1, 0, 0, 0,
              0, 0, 4, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 3, 0, 0, 0, 3, 2, 1,
              0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 3, 1, 2, 3),
              ncol = 22, nrow = 4, byrow = TRUE)

K <- ceiling(M/3)
ngen <- 20
eps = 1e-16
max.iter = 512
alleles <- matrix(c(0, 0,
                    1, 1,
                    2, 2,
                    3, 3,
                    0, 1,
                    0, 2,
                    0, 3,
                    1, 2,
                    1, 3,
                    2, 3), ncol = 2, byrow = T)

# functions
# ============================================================

# function to calculate F under full sib inbreeding
# beginning F = 0
# F[i] = 0.25 + 0.5F[i-1] + 0.25F[i-2]
# ============================================================

calcF <- function(ngen) {
  
  F <- numeric(ngen + 2)
  
  for (i in 1:ngen) {
    F[i + 2] <- 0.25 + 0.5*F[i + 1] + 0.25*F[i]
  }
  
  return(F[ngen+2])
  
}

# calculate likelihood given p and e
# the allele frequency and error vectors
# ============================================================

jgil.ll <- function(R, ngen, p0, e0) {
  
  n.ind <- ncol(R)
  v <- rep(0, 22)
  Fn <- calcF(ngen)
  
  for (i in 1:4) {
    
    v[i] <- p0[i]^2*(1 - Fn) + p0[i]*Fn
    
  }
  i = 5
  for (nt1 in 1:3) {
    for (nt2 in (nt1 + 1):4) {
      
      v[i] <- 0.6 * p0[nt1] * p0[nt2] * (1 - Fn)
      v[i+1] <- 0.8 * p0[nt1] * p0[nt2] * (1 - Fn)
      v[i+2] <- 0.6 * p0[nt1] * p0[nt2] * (1 - Fn)
      i = i + 3
      
    }
  }
  
  A <- matrix(0, ncol = 22, nrow = 4)
  B <- matrix(0, ncol = 22, nrow = 4)
  total.probs <- matrix(0, ncol = 22, nrow = 4)
  S <- matrix(0, ncol = 22, nrow = 4)
  
  err.vec <- 1 - sum(e0) + e0
  
  for (i in 1:22) {
    
    A[, i] <- (4  - M[, i])/4 * e0
    B[, i] <- M[, i]/4 * err.vec
    total.probs[, i] <- A[, i] + B[, i]
    
  }
  
  total.ll <- 0
  
  for (i in 1:n.ind) {
    
    v.idx <- which(v > 0)
    val <- rep(NA, length(v.idx))
    for (j in v.idx) {
      val[j] <- log(v[j]) + sum((log(total.probs[, j])*R[, i])[total.probs[, j] > 0])
    }
    
    mval <- max(val, na.rm = TRUE)
    this.ll <- mval + log(sum(exp(val - mval), na.rm = TRUE))
    total.ll <- total.ll + this.ll
    
  }
  
  return(total.ll)
  
}

# function to maximize likelihood and obtain calls
# ============================================================

jgil <- function(R, ngen, eps = 1e-16, max.iter, M, K, alleles) {
  
  # initialize
  e0 <- rep(0, 4)
  p0 <- rep(0, 4)
  n.iter <- 0
  e1 <- rep(0.05, 4)
  R.sum <- sum(R)
  p1 <- rowSums(R)/R.sum
  n.ind <- ncol(R)
  Fn <- calcF(ngen)
  
  while (sum((p1 - p0)^2) + sum((e1 - e0)^2) > eps & n.iter < max.iter) {
    
    p0 <- p1; e0 <- e1; # update
    
    # genotype frequencies in the vector v based on p0
    v <- rep(0, 22)
    
    for (i in 1:4) {
      
      v[i] <- p0[i]^2*(1 - Fn) + p0[i]*Fn
      
    }
    i = 5
    for (nt1 in 1:3) {
      for (nt2 in (nt1 + 1):4) {
        
        v[i] <- 0.6 * p0[nt1] * p0[nt2] * (1 - Fn)
        v[i+1] <- 0.8 * p0[nt1] * p0[nt2] * (1 - Fn)
        v[i+2] <- 0.6 * p0[nt1] * p0[nt2] * (1 - Fn)
        i = i + 3
        
      }
    }
    
    # calculate e1 and p1, new allele frequency and error
    A <- matrix(0, ncol = 22, nrow = 4)
    B <- matrix(0, ncol = 22, nrow = 4)
    total.probs <- matrix(0, ncol = 22, nrow = 4)
    S <- matrix(0, ncol = 22, nrow = 4)
    
    err.vec <- 1 - sum(e0) + e0
    
    for (i in 1:22) {
      
      A[, i] <- (4  - M[, i])/4 * e0
      B[, i] <- M[, i]/4 * err.vec
      total.probs[, i] <- A[, i] + B[, i]
      
    }
    
    S[total.probs > 0] <- A[total.probs > 0]/total.probs[total.probs > 0]
    
    tmp <- t(R) %*% log(total.probs + eps) + matrix(1, nrow = n.ind) %*% log(v + eps)
    mtmp <- apply(tmp, 1, max)
    J <- exp(tmp - matrix(mtmp, ncol = 1) %*% matrix(1, ncol = 22))
    H <- J/(rowSums(J) %*% matrix(1, ncol = 22))
    
    e1 <- rowSums(S * (R %*% H))/R.sum
    p1 <- K %*% t(H) %*% matrix(1, nrow = n.ind)
    p1 <- p1/sum(p1)
    
    n.iter <- n.iter + 1
    
  }
  
  post.probs <- matrix(0, ncol = 10, nrow = n.ind)
  post.probs[, 1:4] <- H[, 1:4]
  post.probs[, 5] <- rowSums(H[, 5:7])
  post.probs[, 6] <- rowSums(H[, 8:10])
  post.probs[, 7] <- rowSums(H[, 11:13])
  post.probs[, 8] <- rowSums(H[, 14:16])
  post.probs[, 9] <- rowSums(H[, 17:19])
  post.probs[, 10] <- rowSums(H[, 20:22])
  post.probs <- ifelse(post.probs > 1, 1, post.probs) # take care of some cases when the post is a small eps larger than 1
  
  # compute full log likelihood
  full.ll <- jgil.ll(R = R, ngen = ngen, p0 = p1, e0 = e1)
  
  # compute null likelihood, this null is that the lines don't vary
  max.p.i <- which.max(p1)
  p.null <- rep(0, 4)
  p.null[max.p.i] <- 1
  # original jgil method
  # e.null <- rowSums(R)/R.sum
  # e.null[max.p.i] <-0
  e.null <- e1 # however, I think the null should be conditional on the e1, ML estimated error profile
  null.ll1 <- jgil.ll(R = R, ngen = ngen, p0 = p.null, e0 = e.null)
  
  # there is a different null, which is nothing differs from the reference
  p.null <- c(1, 0, 0, 0)
  e.null <- e1
  null.ll2 <- jgil.ll(R = R, ngen = ngen, p0 = p.null, e0 = e.null)
  
  # compute variant quality
  snpQ1 <- min(c(999, floor(-10*log10(pchisq(-2*(null.ll1 - full.ll), 3, lower.tail = F) + 10^-100))))
  snpQ2 <- min(c(999, floor(-10*log10(pchisq(-2*(null.ll2 - full.ll), 3, lower.tail = F) + 10^-100))))
  
  # fetch variant call and call Q
  call.idx <- apply(post.probs, 1, which.max)
  call.alleles <- alleles[call.idx, ]
  callQ <- floor(-10*log10(1 - post.probs[cbind(1:n.ind, call.idx)] + 10^-20))
  callQ <- ifelse(callQ > 99, 99, callQ)
  
  return(list(call = call.alleles, snpQ1 = snpQ1, snpQ2 = snpQ2, callQ = callQ))
  
}

# read lines
# ============================================================

input.con <- file(input, open = "r")
line <- readLines(input.con, n = 1)
ids <- unlist(strsplit(line, split = " "))
n.strain <- length(ids) - 4
write("##fileformat=VCFv4.2", file = output, append = FALSE)
write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file = output, append = TRUE)
write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed, obtained from GATK\">", file = output, append = TRUE)
write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">", file = output, append = TRUE)
write("##INFO=<ID=MDP,Number=1,Type=Float,Description=\"Mean DP\">", file = output, append = TRUE)
write("##INFO=<ID=SDP,Number=1,Type=Float,Description=\"Standard Deviation DP\">", file = output, append = TRUE)
write("##INFO=<ID=JQ,Number=1,Type=Integer,Description=\"Original JGIL snpQ, test for polymorphism\">", file = output, append = TRUE)
write(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", ids[-(1:4)]), file = output, ncolumns =  n.strain + 9, append = TRUE, sep = "\t")

while (TRUE) {
  
  line <- readLines(input.con, n = 1)
  if (length(line) == 0) { break }
  
  allele.count <- unlist(strsplit(line, split = " "))
  snp.info <- allele.count[1:4]
  r.count <- matrix(as.integer(unlist(strsplit(allele.count[5:(n.strain + 4)], split = ","))), ncol = n.strain)
  r.mean <- round(mean(colSums(r.count)), 2)
  r.sd <- round(sd(colSums(r.count)), 2)
  R.count <- matrix(0, nrow = 4, ncol = n.strain)
  R.count[1:nrow(r.count), ] <- r.count
  this.jgil <- jgil(R = R.count, ngen = 20, eps = 1e-16, max.iter = 512, M = M, K = K, alleles = alleles)
  write(c(snp.info[1:2], ".", snp.info[3:4], this.jgil$snpQ2, ".",
          paste("MDP=", r.mean, ";SDP=", r.sd, ";JQ=", this.jgil$snpQ1, sep = ""),
          "GT:AD:GQ", paste(apply(this.jgil$call, 1, paste, collapse = "/"),
          apply(r.count, 2, paste, collapse = ","), this.jgil$callQ, sep = ":")),
          file = output, ncolumns = n.strain + 9, append = TRUE, sep = "\t")

}

close(input.con)
