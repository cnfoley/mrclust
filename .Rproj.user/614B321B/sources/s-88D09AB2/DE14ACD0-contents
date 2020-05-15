clust_prob <- function(bic_clust, clust_mn, clust_pi, theta, theta_sd, df,
                       junk_mixture = FALSE, junk_mean = NULL, junk_sd = NULL,
                       null_mixture = FALSE, null_mean = NULL, null_sd = NULL,
                       obs_names, clust_size_prior = FALSE, prior = NULL) {
  if (!clust_size_prior) {
    mx <- which.min(bic_clust)
  } else {
    m <- length(theta)
    pr <- stats::dbinom(1:m, m, prior)
    pr <- pr[seq_len(bic_clust)] / sum(pr[seq_len(bic_clust)])
    mx <- which.max(exp(-bic_clust / 2) * pr)
  }
  theta_mn <- clust_mn[[mx]]
  clust_probs <- clust_pi[[mx]]
  n_clust <- length(theta_mn)
  m <- length(theta)
  res1 <- round(as.numeric(t(rij(1:m, 1:n_clust, clust_probs,
                                     theta, theta_sd, theta_mn,
                                     junk_mixture, df, junk_mean, junk_sd,
                                     null_mixture, null_mean, null_sd))), 3)
  res2 <- round(rep(theta_mn, m), 3)
  res <- data.table::data.table(rep(obs_names, each = n_clust),
                                rep(1:n_clust, m),
                                res1,
                                res2,
                                rep(theta, each = n_clust),
                                rep(theta_sd, each = n_clust))
  names(res) <- c("observation", "cluster", "probability",
                  "cluster_mean", "theta", "theta_se")
  factor_clust <- res$cluster
  if (null_mixture & !junk_mixture) {
    tmp <- res$cluster == n_clust
    factor_clust[tmp] <- "Null"
  } else if (null_mixture & junk_mixture) {
    tmp <- res$cluster == n_clust - 1
    factor_clust[tmp] <- "Null"
  }
  if (junk_mixture) {
    tmp <- res$cluster == n_clust
    factor_clust[tmp] <- "Junk"
  }
  res$cluster_class <- as.character(factor_clust)
  return(res)
}
