rij <- function(i, j, pi_clust, theta, theta_sd, theta_clust, junk_mixture,
                df, mu, sig, null_mixture, mu_null, sig_null) {
  tmp_lgt <- length(i)
  if (junk_mixture | null_mixture) {
    tmp_lgt2 <- length(j) - sum(junk_mixture) - sum(null_mixture)
    tmp_mn <- theta_clust[1:tmp_lgt2]
    j <- 1:(length(j) - sum(junk_mixture) - sum(null_mixture))
  } else {
    tmp_lgt2 <- length(j)
    tmp_mn <- theta_clust
  }
  if (junk_mixture | null_mixture) {
    tmp_lg_num <- log_norm(i, j, tmp_lgt, tmp_lgt2, theta, theta_sd, tmp_mn)
    if (null_mixture) {
      tmp_lg_num <- c(tmp_lg_num, null_den(x = theta, mu = mu_null,
                                           sig = sig_null, log = TRUE))
    }
    if (junk_mixture) {
      tmp_lg_num <- c(tmp_lg_num,
                      gen_t(theta, df = df, mu = mu, sig = sig, log = TRUE))
    }
    num <- exp(matrix(rep(log(pi_clust), each = tmp_lgt) + tmp_lg_num,
                      nrow = tmp_lgt,
                      ncol = tmp_lgt2 + sum(junk_mixture) + sum(null_mixture)))
    num[num < 1e-300] <- 1e-300
    den <- rowSums(num)
  } else {
    tmp_lg_num <- log_norm(i, j = seq_len(length(theta_clust)), tmp_lgt,
                           tmp_lgt2 = length(theta_clust),
                           theta, theta_sd, theta_clust)
    num <- exp(matrix(rep(log(pi_clust[j]), each = tmp_lgt) + tmp_lg_num,
                      nrow = tmp_lgt, ncol = tmp_lgt2))
    num[num < 1e-300] <- 1e-300
    den <- rowSums(num)
  }
  res <- num / den
  return(res)
}
