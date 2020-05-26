pi_updt <- function(j, m, pi_clust, theta, theta_sd, theta_clust, junk_mixture,
                    df, mu, sig, junk_prob, fix_junk_prob, null_mixture,
                    mu_null, sig_null, null_prob, fix_null_prob) {
  tmp <- rij(1:m, j, pi_clust, theta, theta_sd, theta_clust, junk_mixture,
             df, mu, sig, null_mixture, mu_null, sig_null)
  num <- if (length(j) > 1) {
    colSums(tmp)
  } else {
    sum(tmp)
  }
  jk_null <- sum(junk_mixture) + sum(null_mixture)
  n_prob <- j_prob <- FALSE
  if (null_mixture & !is.null(null_prob)) {
    tmp_lgt <- length(num) - sum(junk_mixture)
    if (num[tmp_lgt] / m < null_prob | fix_null_prob == TRUE) {
      num[tmp_lgt] <- null_prob * m
      n_prob <- TRUE
    }
  }
  if (junk_mixture & !is.null(junk_prob)) {
    tmp_lgt <- length(num)
    if (num[length(num)] / m < junk_prob | fix_junk_prob == TRUE) {
      num[tmp_lgt] <- junk_prob * m
      j_prob <- TRUE
    }
  }
  if (j_prob | n_prob) {
    num[1:(length(num) - jk_null)] <- (m *
      (1 - sum(num[(length(num) - jk_null + 1):length(num)]) / m) *
      num[1:(length(num) - jk_null)] / (sum(num[1:(length(num) - jk_null)])))
  }
  return(num / m)
}
