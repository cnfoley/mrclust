theta_updt <- function(j, m, theta, theta_sd, theta_clust, pi_clust,
                       junk_mixture, df, mu, sig, null_mixture, mu_null,
                       sig_null) {
  tmp <- rij(1:m, j, pi_clust, theta, theta_sd, theta_clust, junk_mixture, df,
             mu, sig, null_mixture, mu_null, sig_null) / theta_sd[1:m]^2
  den <- if (length(j) > 1) {
    colSums(tmp)
  } else {
    sum(tmp)
  }
  num <- if (length(j) > 1) {
    colSums(tmp * theta)
  } else {
    sum(tmp * theta)
  }
  new_theta <- num / den
  null_jk <- sum(null_mixture) + sum(junk_mixture)
  if (null_mixture & null_jk == 2) {
    new_theta[length(j) - 1] <- mu_null
  } else if (null_mixture & null_jk == 1) {
    new_theta[length(j)] <- mu_null
  }
  if (junk_mixture) {
    new_theta[length(j)] <- mu
  }
  return(num / den)
}
