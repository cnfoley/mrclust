#' Cluster size and assignemnt probabilities
#'
#' Keep results based on a minimum allocation probability and number of observations in a cluster.
#'
#' @param dta table of results from mr_clust_em$results$best.
#' @param prob numeric scalar, keep only variants assigned to clusters above this allocation probability.
#' @param min.obs integer, keep only variants assinged to clusters with more than or equal to min.obs members.
#' @return The results
#' @export

pr.clust <- function(dta, prob = 0.5, min.obs = 1) {
  dta <- dta[dta$probability > prob, ]
  tmp <- unique(dta$cluster_mean)
  tmp2 <- rep(1, dim(dta)[1])
  for (i in tmp) {
    tmp3 <- (dta$cluster_mean == i)
    if (sum(tmp3) > 1) {
      tmp2[tmp3] <- sum(tmp3)
    }
  }
  dta$cluster.size <- tmp2
  return(dta[dta$cluster.size >= min.obs, ])
}
