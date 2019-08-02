# D.L. Giraldo (2019)

# Calculate Haralick features from a Gray-level Co-ocurrece matrix (GLCM) 
# these calculations assume the GLCM is symmetric
Haralick_14features <- function(G){
  p <- G/sum(G)
  
  ### Pre-calculations ###
  # Marginal probability (p_x = p_y because symmetry)
  i <- seq(1, nrow(p))
  p_m <- rowSums(p)
  mu_m <- sum(i*p_m)
  var_m <- sum((i-mu_m)^2*p_m)
  # Difference probability
  k <- seq(0, nrow(p)-1)
  p_diff <- sapply(k, function(x) sum(p[abs(row(p) - col(p)) == x]))
  # Sum probability
  l <- seq(2, 2*nrow(p))
  p_sum <- sapply(l, function(x) sum(p[(row(p) + col(p)) == x]))
  # Entropy
  HXY <- -sum(p[p > 0]*log2(p[p > 0]))
  
  HF <- vector("double", length = 14)
  # ASM = Energy
  HF[1] <- sum(p^2)
  # Contrast
  HF[2] <- sum(k^2*p_diff)
  # Correlation
  # HF[3] <- sum(row(p)*col(p)*p - mu_m^2)/(var_m) <- As it is in Haralick paper.
  HF[3] <- sum((row(p) - mu_m)*(col(p) - mu_m)*p)/var_m
  # Sum of Squares = Variance of marginal distribution
  HF[4] <- var_m
  # Inverse different moment
  HF[5] <- sum(p/(1+(row(p)-col(p))^2))
  # Sum Average
  HF[6] <- sum(l*p_sum)
  # Sum Variance
  HF[7] <- sum((l-sum(l*p_sum))^2*p_sum)
  # Sum Entropy
  HF[8] <- -sum(p_sum[p_sum > 0]*log2(p_sum[p_sum > 0]))
  # Entropy
  HF[9] <- HXY
  # Difference Variance
  HF[10] <- sum((k-sum(k*p_diff))^2*p_diff)
  # Difference Entropy
  HF[11] <- -sum(p_diff[p_diff > 0]*log2(p_diff[p_diff > 0])) 
  # Information Measures of Correlation
  HXY1 <- -sum(p*log2((p_m %o% p_m) + 1e-6)) 
  H_m <- -sum(p_m[p_m > 0]*log2(p_m[p_m > 0]))
  HF[12] <- (HXY - HXY1)/H_m
  HXY2 <- -sum((p_m %o% p_m)*log2((p_m %o% p_m) + 1e-6))
  HF[13] <- sqrt(1 - exp(-2*(HXY2 - HXY)))
  # Maximal Correlation Coeficient
  Q <- array(0, dim = dim(p))
  for (n in i){
    for (m in i){
      Q[n,m] <- sum((p[n,i]*p[m,i])/(p_m[n]*p_m), na.rm = TRUE)
    }
  }
  HF[14] <- sqrt(eigen(Q, symmetric = FALSE, only.values = TRUE)$values[2])
  return(HF)
}