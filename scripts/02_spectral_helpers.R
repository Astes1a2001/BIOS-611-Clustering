
## ---------- helper: adjacency ----------
build_adjacency <- function(x, d_threshold, self_loops = FALSE) {
  X <- as.matrix(x)
  D <- as.matrix(dist(X))
  A <- (D < d_threshold) * 1L
  diag(A) <- if (self_loops) 1L else 0L
  storage.mode(A) <- "double"
  A
}

## ---------- helper: symmetric normalized Laplacian ----------
laplacian_sym <- function(A) {
  n <- nrow(A)
  deg <- rowSums(A)
  Dinvhalf <- diag(ifelse(deg > 0, 1 / sqrt(deg), 0))
  Lsym <- diag(n) - Dinvhalf %*% A %*% Dinvhalf
  Lsym
}

## ---------- helper: eigen-decomposition & pick k smallest ----------
laplacian_eigen <- function(Lsym, k) {
  eg <- eigen(Lsym, symmetric = TRUE)
  ord <- order(eg$values, decreasing = FALSE)
  list(
    U = eg$vectors[, ord[seq_len(k)], drop = FALSE],
    values = eg$values[ord[seq_len(k)]]
  )
}

## ---------- helper: row-normalize & kmeans ----------
cluster_eigenvectors <- function(U, k, nstart = 20, iter.max = 50, row_normalize = TRUE) {
  Y <- U
  if (row_normalize) {
    rnorms <- sqrt(rowSums(U^2))
    rnorms[rnorms == 0] <- 1
    Y <- U / rnorms
  }
  km <- kmeans(Y, centers = k, nstart = nstart, iter.max = iter.max)
  km$cluster
}

## ---------- spectral clustering wrapper compatible with clusGap ----------
spectral_kmeans <- function(x,
                            k,
                            d_threshold = 1.0,
                            self_loops = FALSE,
                            nstart = 20,
                            iter.max = 50,
                            row_normalize = TRUE) {
  A    <- build_adjacency(x, d_threshold = d_threshold, self_loops = self_loops)
  Lsym <- laplacian_sym(A)
  eig  <- laplacian_eigen(Lsym, k = k)
  U    <- eig$U
  
  labels <- cluster_eigenvectors(
    U,
    k = k,
    nstart = nstart,
    iter.max = iter.max,
    row_normalize = row_normalize
  )
  
  list(cluster = labels)
}
