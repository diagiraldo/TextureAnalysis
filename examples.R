source("~/MyTextureAnalysis/glcm_dir.R")
source("~/MyTextureAnalysis/Haralick_14features.R")

# Simulate random 2D images
set.seed(1987)
im_dim <- 2
l <- 10
a1 <- array(runif(l^im_dim), dim = rep(l, im_dim))
a2 <- row(a1)*col(a1)
a2 <- a2[, ncol(a2):1]
a3 <- 0*a2
a3[lower.tri(a3, diag = TRUE)] <- 1
a3[lower.tri(a3)] <- 2
a4 <- a3 + a1

A <- list(a1, a2, a3, a4)

# Plot those images
for (i in 1:length(A)){
  png(sprintf("~/MyTextureAnalysis/ex%d_image2.png", i))
  par(mar = c(0,0,0,0))
  image(t(A[[i]])[,nrow(A[[i]]):1], axes = FALSE, col = gray.colors(32))
  dev.off()
}

# List of quantized images
n_q <- 4
Q <- list()
for (i in 1:length(A)){
  br <- quantile(as.vector(A[[i]]), probs = seq(0,1,1/n_q), names = FALSE)
  I <- array(cut(A[[i]], unique(br), labels = FALSE, include.lowest = TRUE), dim = dim(A[[i]]))
  Q[[i]] <- I
  png(sprintf("~/MyTextureAnalysis/ex%d_quant2.png", i))
  par(mar = c(0,0,0,0))
  image(t(I)[,nrow(I):1], axes = FALSE, col = gray.colors(n_q))
  dev.off()
}

# Parameters for GLCM
distance <- 1
dirlist <- list(c(0, 1), c(-1, 1), c(-1, 0), c(-1, -1))
n_dirs = length(dirlist);

# List of GLCM
G <- list()
for (i in 1:length(Q)){
  glcm = array(0, dim = c(n_q, n_q, n_dirs))
  for (d in 1:n_dirs){
    glcm[,,d] <- glcm_dir(Q[[i]], n_q, distance, dirlist[[d]], symmetric = TRUE)
  }
  G[[i]] <- apply(glcm, c(1,2), sum)
}

# Haralick Features
harnames <- c("Angular Second Moment", "Contrast", "Correlation", "Variance", "Inverse Difference Moment",
              "Sum Average", "Sum Variance", "Sum Entropy", "Entropy", "Difference Variance", "Difference Entropy",
              "Info. Measure of Corr. I", "Info. Measure of Corr. II", "Maximal Corr. Coeff.")
H <- data.frame(row.names = harnames)
for (i in 1:length(G)){
  H[,i] <- Haralick_14features(G[[i]])
}
H <- round(H, digits = 2)
H <- mutate(H, feat = harnames) %>%
  dplyr::select(feat, V1, V2, V3, V4)
paste(apply(H, 1, function(x) paste(x, collapse = " & ")), collapse = "\\ ")

# Angular Haralick
i <- 4
G <- list()
H <- data.frame(row.names = harnames)
for (d in 1:n_dirs){
  G[[d]] <- glcm_dir(Q[[i]], n_q, distance, dirlist[[d]], symmetric = TRUE)
  H[,d] <- Haralick_14features(G[[d]])
}
H <- round(H, digits = 2)
H <- mutate(H, feat = harnames) %>%
  dplyr::select(feat, V1, V2, V3, V4) %>%
  mutate(mean = round((V1+V2+V3+V4)/4, digits = 2), range = pmax(V1, V2, V3, V4)- pmin(V1, V2, V3, V4))
paste(apply(H, 1, function(x) paste(x, collapse = " & ")), collapse = "\\ ")


##### Print matrix as latex
m <- cbind(1:n_q, G[[2]])
paste(paste(" ", 1:n_q, collapse = " & "), paste(apply(m, 1, function(x) paste(x, collapse = " & ")), collapse = "\\ "), collapse = "\\")

g <- glcm_dir(Q[[1]], 4, 1, c(0,1), symmetric = FALSE) 
paste(apply(g, 1, function(x) paste(x, collapse = " & ")), collapse = "\\ ")

library(dplyr)
H <- read.table("~/MyTextureAnalysis/HF4presentation.txt", stringsAsFactors = FALSE, header = FALSE, sep = ";")
H <- round(H, digits = 2)

H <- mutate(H, feat = harnames) %>%
  dplyr::select(feat, V1, V2, V3, V4)

