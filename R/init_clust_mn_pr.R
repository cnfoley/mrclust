init_clust_mn_pr <- function(theta, range = length(theta) - 1, iters = 1e4) {
  dta <- vector("list", iters)
  dta_mn <- vector("list", iters)
  dta_bic <- vector("numeric", iters)
  for (i in 1:iters) {
    tmp_dta <- theta
    if (length(range) > 1) {
      num_c <- sample(range, 1)
    } else {
      num_c <- range
    }
    tmp_res <- stats::kmeans(tmp_dta, num_c)
    dta[[i]] <- tmp_res$cluster
    dta_mn[[i]] <- as.numeric(tmp_res$centers)
    dta_bic [i] <- kmeansBIC(tmp_res)
  }
  mx <- which.min(dta_bic)
  num_c <- length(unique(dta[[mx]]))
  tmp_mn <- as.numeric(dta_mn[[mx]])
  tmp_clust <- dta[[mx]]

  m <- length(theta)
  tmp_dta <- sample(theta, m)
  tmp_res <- stats::kmeans(tmp_dta, num_c)
  tmp_pr <- table(tmp_res$cluster)

  for (i in 1:(iters - 1)) {
    tmp_dta <- sample(theta, m)
    tmp_res <- stats::kmeans(tmp_dta, num_c)
    tmp_pr <- tmp_pr + table(tmp_res$cluster)
  }
  tmp_pr <- tmp_pr / (m * iters)

  return(list(tmp_clust, tmp_mn, tmp_pr))
}
