clust_inc_list <- function(dta, by_prob = 0.2, bound = 0) {
  if (!((1 - bound) / by_prob) %% 1 == 0) {
    stop("(1-bound)/by_prob needs to be an integer")
  }
  unique_clust <- unique(dta$cluster)
  res <- vector("list", length(unique_clust))
  clust_nam <- vector("character", length(unique_clust))
  nms_tmp <- unique(dta$cluster_class)
  for (i in unique_clust) {
    if (bound == 0) {
      sub_grps <- (1 / by_prob)
    } else {
      sub_grps <- ceiling((1 - bound) / by_prob)
    }
    tmp <- vector("list", sub_grps)
    tmp_nam <- vector("character", sub_grps)
    for (j in 1:sub_grps) {
      p_mx <- 1 - (j - 1) * by_prob
      p_mn <- 1 - j * by_prob
      if (p.mx <= bound) {
        break
      }
      if ((!bound == 0) & (1 - j * by_prob < bound)) {
        p_mn <- bound
      }
      tmp[[j]] <- (dta$observation[(dta$probability <= p_mx)
                                   & (dta$probability > p_mn)
                                   & (dta$cluster == i)])
      tmp_nam[j] <- paste0("pr_range_", p_mx, "_to_", p_mn)
    }
    names(tmp) <- tmp_nam
    clust_nam[i] <- paste0("cluster_", nms_tmp[i])
    res[[i]] <- tmp
  }
  names(res) <- clust_nam
  return(res)
}
