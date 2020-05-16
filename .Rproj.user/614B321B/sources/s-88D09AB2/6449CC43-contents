#' Plotting clustered ratio-estimates
#'
#' Plot of the two-stage regression estimates, i.e. G-X and G-Y associations,
#' annotated with cluster allocation labels and cluster mean estimates.
#'
#' @param res table of results from mr_clust_em$results$best.
#' @param bx numeric vector of length the number of variants, the i-th element
#' is the estimated regression coefficient - i.e. beta-x value - relating the
#' i-th genetic variant to the risk-factor.
#' @param by numeric vector of length the number of variants, the i-th element
#' is the estimated regression coefficient - i.e. beta-y value - relating the
#' -th genetic variant to the outcome.
#' @param bxse,byse numeric vector of length the number of variants, the i-th
#' element is the standard error of the estimated regression coefficient
#' relating the i-th genetic variant to the risk-factor.
#' @param obs_names character vector of length the number of variants, the i-th
#' element is the name of the i-th genetic variants - e.g. the rsID.
#' @return Returned is a scatter plot of the two-stage association estimates
#' for each variant in which: clusters are colour coded and variants with larger
#'  assignement/inclusion probabilities appear larger.
#' @export

two_stage_plot <- function(res, bx, by, bxse, byse, obs_names) {
  cbpalette <- c("Null" = "#CC79A7", "Junk" = "#000000", "1" = "#999999",
                 "2" = "#0072B2", "3" = "#D55E00", "4" = "#F0E442",
                 "5" = "#009E73", "6" = "#56B4E9", "7" = "#E69F00")
  clst_class <- unique(res$cluster_class)
  no_clst <- clst_class[!clst_class %in% c("Null", "Junk")]
  if (!identical(no_clst, character(0))) {
    n_palette <- max(as.numeric(no_clst))
  } else {
    n_palette <- 0
  }
  if (n_palette >= 8) {
    tmp_cbpalette <- RColorBrewer::brewer.pal(max(3, n_palette),
                                              "Spectral")[8:n_palette]
    names(tmp_cbpalette) <- paste0(8:n_palette)
    cbpalette <- c(cbpalette, tmp_cbpalette)
  }

  tmp_orient <- orientate(bx, by)
  bx <- tmp_orient$bx
  by <- tmp_orient$by

  # redefine null cluster as cluster zero
  res$cluster[res$cluster_mean == 0] <- 0
  res$clusters <- as.factor(res$cluster)

  ord <- match(res$observation, obs_names)
  res$bx <- bx[ord]
  res$by <- by[ord]
  res$bxse <- bxse[ord]
  res$byse <- byse[ord]

  p <- ggplot2::ggplot(data = res, ggplot2::aes(bx, by)) +
    ggplot2::geom_point(ggplot2::aes(colour = cluster_class,
                                     size = probability), shape = 21) +
    ggplot2::scale_color_manual(values = cbpalette) +
    ggplot2::labs(size = "Cluster\ninclusion\nprobability") +
    ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = bx - 1.96 * bxse, xmax = bx + 1.96 * bxse,
                   color = cluster_class), linetype = "solid") +
    ggplot2::geom_errorbar(
    ggplot2::aes(ymin = by - 1.96 * byse, ymax = by + 1.96 * byse,
                   color = cluster_class), linetype = "solid") +
    ggplot2::ylab("Genetic effect on outcome") +
    ggplot2::xlab("Genetic effect on exposure") +
    ggplot2::ggtitle("Clustered two-stage ratio estimates") +
    ggplot2::theme_bw()

  clusts <- unique(res$cluster_class)
  clusts <- clusts[order(clusts)]
  ind <- which(!clusts %in% "Null" & !clusts %in% "Junk")
  if (sum(ind) > 0) {
    clust_label <- c("Junk cluster")
  } else {
    clust_label <- NULL
  }
  if (sum(ind) > 0) {
    for (i in seq_len(ind)) {
      mn <- res$cluster_mean[res$cluster == as.numeric(clusts[ind[i]])][1]
      clust_label <- c(clust_label, paste0("cluster_", ind[i]))
      col <- cbpalette[names(cbpalette) == clusts[ind[i]]]
      p <- p + ggplot2::geom_abline(slope = mn, intercept = 0, color = col,
                                    linetype = "dotted", lwd = 0.50)
    }
  }
  p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Cluster",
                                                         labels = clust_label))
  p + ggplot2::xlim(0, 1.1 * max(abs(bx)))
}
