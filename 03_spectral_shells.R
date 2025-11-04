library(cluster)
library(ggplot2)

source("scripts/02_spectral_helpers.R")

generate_shell_clusters <- function(n_shells, k_per_shell, max_radius,
                                    noise_sd = 0.1, min_radius_fraction = 0.2) {
  min_radius <- min_radius_fraction * max_radius
  radii_raw <- ((1:n_shells)^2)
  radii_raw <- (radii_raw - min(radii_raw)) / (max(radii_raw) - min(radii_raw))
  radii_raw <- min_radius + radii_raw * (max_radius - min_radius)
  
  total_points <- n_shells * k_per_shell
  xs <- ys <- zs <- numeric(total_points)
  labels <- integer(total_points)
  idx_start <- 1
  
  for (shell_id in seq_len(n_shells)) {
    base_r <- radii_raw[shell_id]
    r_vec <- pmax(base_r + rnorm(k_per_shell, 0, noise_sd), 1e-6)
    theta <- runif(k_per_shell, 0, 2*pi)
    u <- runif(k_per_shell, -1, 1)
    phi <- acos(u)
    idx_end <- idx_start + k_per_shell - 1
    xs[idx_start:idx_end] <- r_vec * sin(phi) * cos(theta)
    ys[idx_start:idx_end] <- r_vec * sin(phi) * sin(theta)
    zs[idx_start:idx_end] <- r_vec * cos(phi)
    labels[idx_start:idx_end] <- shell_id
    idx_start <- idx_end + 1
  }
  
  data.frame(x = xs, y = ys, z = zs, shell = labels)
}

set.seed(42)

n_shells    <- 4
k_per_shell <- 100
noise_sd    <- 0.1
d_threshold <- 1

radii_to_test <- seq(10, 0, by = -1)

results <- data.frame(
  max_radius = numeric(0),
  k_hat      = integer(0)
)

for (R in radii_to_test) {
  dat <- generate_shell_clusters(
    n_shells    = n_shells,
    k_per_shell = k_per_shell,
    max_radius  = R,
    noise_sd    = noise_sd
  )
  
  X <- as.matrix(dat[, c("x","y","z")])
  
  gap <- clusGap(
    X,
    FUN = spectral_kmeans,
    K.max = n_shells,
    B = 20,
    d_threshold = d_threshold,
    nstart = 20,
    iter.max = 50
  )
  
  k_hat <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
  
  results <- rbind(results, data.frame(
    max_radius = R,
    k_hat      = k_hat
  ))
}

dir.create("figs", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

png("figs/spectral_shells_gap.png", width = 700, height = 500)
plot(results$max_radius,
     results$k_hat,
     type = "b",
     xlab = "max_radius (outer shell size)",
     ylab = "Estimated # of clusters (k_hat)",
     main = "Spectral clustering recoverability vs shell separation")
abline(h = n_shells, lty = 2)
text(x = max(results$max_radius),
     y = n_shells,
     labels = paste0("true k = ", n_shells),
     pos = 3)
dev.off()

write.csv(results, "results/spectral_shells_results.csv", row.names = FALSE)


