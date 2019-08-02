# D.L. Giraldo (2019)
# Calculate the Gray-level Co-ocurrece matrix (GLCM) given a quantized image 

glcm_dir <- function(I, n_levels, distance, dir, symmetric = TRUE){
  n_q <- n_levels
  im_dim <- length(dim(I))
  G = array(0, dim = c(n_q, n_q))
  for (q in 1:n_q){
    idx <- which(I == q, arr.ind = TRUE)
    if (nrow(idx) > 0){
      idx_dir <- t(apply(idx, 1, function(x) x+distance*dir))
      idx_dir <- subset(idx_dir, (apply(idx_dir, 1, min) > 0 & apply(idx_dir <= dim(I), 1, sum) == im_dim))
      if (nrow(idx_dir) >= 1){
        for (c in 1:nrow(idx_dir)){
          valq <- ifelse(im_dim == 3, I[idx_dir[c,1], idx_dir[c,2], idx_dir[c,3]], I[idx_dir[c,1], idx_dir[c,2]])
          G[q, valq] <- G[q, valq] + 1
          if (symmetric) {G[valq, q] <- G[valq, q] + 1}
        }
      }
    }
  }
  return(G)
}



