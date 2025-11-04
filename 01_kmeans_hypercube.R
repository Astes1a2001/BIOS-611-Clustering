library(cluster)
library(ggplot2)
library(plotly)

generate <- function(n, k, side_length, noise_sd = 1.0) {
  centers <- diag(side_length, nrow = n, ncol = n)
  make_cluster <- function(center) {
    matrix(rep(center, each = k), nrow = k) +
      matrix(rnorm(k * n, 0, noise_sd), nrow = k)
  }
  X <- do.call(rbind, lapply(1:n, function(i) make_cluster(centers[i, ])))
  as.data.frame(X)
}

set.seed(111)
dims <- c(6, 5, 4, 3, 2)
side_lengths <- 10:1
k <- 100
noise_sd <- 1.0

results <- data.frame()

for (n in dims) {
  for (L in side_lengths) {
    X <- generate(n, k, L, noise_sd)
    gap <- clusGap(as.matrix(X),
                   FUN = kmeans,
                   K.max = n,
                   nstart = 20,
                   iter.max = 50)
    k_hat <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method = "firstSEmax")
    results <- rbind(results, data.frame(n = n, side_length = L, k_est = k_hat))
  }
}

p <- ggplot(results, aes(x = side_length, y = k_est, color = factor(n))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = n), linetype = "dashed", color = "gray40") +
  scale_x_reverse(breaks = 10:1) +
  labs(
    title = "Estimated Number of Clusters vs. Side Length",
    x = "Side Length (L)",
    y = "Estimated Clusters (Gap Statistic)",
    color = "Dimension (n)"
  ) +
  theme_minimal(base_size = 13)

dir.create("figs", showWarnings = FALSE)
ggsave("figs/kmeans_hypercube_gap.png", p, width = 7, height = 5, dpi = 300)

write.csv(results, "results/kmeans_hypercube_results.csv", row.names = FALSE)
