loglik <- function(m, pi_clust, theta, theta_sd, theta_clust,
                   junk_mixture, df, mu, sig,
                   null_mixture, mu_null, sig_null) {
  i <- 1:m
  if (junk_mixture | null_mixture) {
    k <- length(theta_clust) - sum(junk_mixture) - sum(null_mixture)
  } else {
    k <- length(theta_clust)
  }
  tmp_lg_den <- log_norm(i, j = 1:k, m, tmp_lgt2 = k, theta, theta_sd,
                         theta_clust)
  if (null_mixture) {
    tmp_lg_den <- c(tmp_lg_den, null_den(x = theta, mu = mu_null,
                                          sig = sig_null, log = TRUE))
  }
  if (junk_mixture) {
    tmp_lg_den <- c(tmp_lg_den, gen_t(theta, df = df, mu = 0, sig = sig,
                                       log = TRUE))
  }

  tmp <- exp(matrix(rep(log(pi_clust), each = m) + tmp_lg_den, nrow = m,
                    ncol = length(pi_clust)))
  tmp_rowsum <- rowSums(tmp)
  res <- sum(log(tmp_rowsum))
  return(res)
}
